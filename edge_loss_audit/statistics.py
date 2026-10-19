"""Core metric and inference functions extracted from the study analysis code.

The calculation bodies are unchanged; source paths are in docs/provenance.md.
"""
from __future__ import annotations
from typing import Any
import numpy as np
import pandas as pd

def per_group_statistics(frame: pd.DataFrame, score_column: str, group_columns: list[str]) -> pd.DataFrame:
    rows = []
    grouper: str | list[str] = group_columns[0] if len(group_columns) == 1 else group_columns
    for key, group in frame.groupby(grouper, sort=False, dropna=False):
        target = group["analysis_target"].to_numpy(dtype=int)
        if np.unique(target).size != 2:
            continue
        score = group[score_column].to_numpy(dtype=float)
        difference = score[target == 1][:, None] - score[target == 0][None, :]
        correct = float((difference > 0).sum() + 0.5 * (difference == 0).sum())
        pairs = int(difference.size)
        rows.append(
            {
                "group": "|".join(map(str, key)) if isinstance(key, tuple) else str(key),
                "macro_auroc": correct / pairs,
                "pair_correct": correct,
                "pairs": pairs,
                "rows": int(len(group)),
            }
        )
    return pd.DataFrame(rows)

def adjust_bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted_ranked = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return adjusted

def adjust_holm(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted_ranked = np.maximum.accumulate(
        ranked * np.arange(len(values), 0, -1)
    )
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return adjusted

def paired_inference(
    first: pd.DataFrame,
    second: pd.DataFrame,
    metric: str,
    n_bootstrap: int,
    n_permutations: int,
    seed: int,
) -> dict[str, Any]:
    merged = first.merge(
        second,
        on="group",
        suffixes=("_first", "_second"),
        validate="one_to_one",
    )
    if not np.array_equal(
        merged["pairs_first"].to_numpy(), merged["pairs_second"].to_numpy()
    ):
        raise RuntimeError("Pair denominators differ between methods")
    if "resampling_cluster_first" in merged:
        if not merged["resampling_cluster_first"].astype(str).equals(
            merged["resampling_cluster_second"].astype(str)
        ):
            raise RuntimeError("Resampling clusters differ between methods")
        merged["resampling_cluster"] = merged["resampling_cluster_first"].astype(str)
    else:
        merged["resampling_cluster"] = merged["group"].astype(str)
    matched_strata = len(merged)
    if metric == "macro_auroc":
        merged["_delta"] = (
            merged["macro_auroc_first"].to_numpy(dtype=float)
            - merged["macro_auroc_second"].to_numpy(dtype=float)
        )
        clustered = merged.groupby("resampling_cluster", sort=False).agg(
            delta=("_delta", "mean")
        )
        delta = clustered["delta"].to_numpy(dtype=float)
        denominator = None
        observed = float(delta.mean())
    elif metric == "pair_accuracy":
        merged["_delta"] = (
            merged["pair_correct_first"].to_numpy(dtype=float)
            - merged["pair_correct_second"].to_numpy(dtype=float)
        )
        clustered = merged.groupby("resampling_cluster", sort=False).agg(
            delta=("_delta", "sum"),
            pairs=("pairs_first", "sum"),
        )
        delta = clustered["delta"].to_numpy(dtype=float)
        denominator = clustered["pairs"].to_numpy(dtype=float)
        observed = float(delta.sum() / denominator.sum())
    else:
        raise ValueError(metric)

    n_groups = len(delta)
    probability = np.full(n_groups, 1.0 / n_groups)
    rng = np.random.default_rng(seed)
    bootstrap = []
    for start in range(0, n_bootstrap, 500):
        size = min(500, n_bootstrap - start)
        counts = rng.multinomial(n_groups, probability, size=size)
        if metric == "macro_auroc":
            bootstrap.append(counts @ delta / n_groups)
        else:
            bootstrap.append((counts @ delta) / (counts @ denominator))
    bootstrap_values = np.concatenate(bootstrap)

    extreme = 0
    completed = 0
    for start in range(0, n_permutations, 1000):
        size = min(1000, n_permutations - start)
        signs = rng.choice(np.array([-1.0, 1.0]), size=(size, n_groups))
        if metric == "macro_auroc":
            null = signs @ delta / n_groups
        else:
            null = (signs @ delta) / denominator.sum()
        extreme += int(np.sum(np.abs(null) >= abs(observed)))
        completed += size
    p_value = float((extreme + 1) / (completed + 1))
    return {
        "matched_strata": matched_strata,
        "groups": n_groups,
        "positive_negative_pairs": int(merged["pairs_first"].sum()),
        "delta": observed,
        "ci_low": float(np.quantile(bootstrap_values, 0.025)),
        "ci_high": float(np.quantile(bootstrap_values, 0.975)),
        "permutation_p": p_value,
        "probability_gt_zero": float(np.mean(bootstrap_values > 0)),
    }

def shuffle_pair_accuracy(frame: pd.DataFrame, score_column: str, n_permutations: int, seed: int) -> dict[str, object]:
    groups = []
    for _, group in frame.groupby("mutation", sort=False):
        if len(group) != 2 or group.target.nunique() != 2:
            continue
        target = group.target.to_numpy(dtype=int)
        score = group[score_column].to_numpy(dtype=float)
        positive_score = float(score[target == 1][0])
        negative_score = float(score[target == 0][0])
        groups.append((positive_score, negative_score))
    if not groups:
        raise ValueError("No two-partner mixed-outcome groups for shuffle control")
    pairs = np.asarray(groups, dtype=float)
    observed = float(np.mean(np.where(pairs[:, 0] > pairs[:, 1], 1.0, np.where(pairs[:, 0] == pairs[:, 1], 0.5, 0.0))))
    rng = np.random.default_rng(seed)
    null = np.empty(n_permutations, dtype=float)
    for index in range(n_permutations):
        swap = rng.integers(0, 2, size=len(pairs)).astype(bool)
        first = np.where(swap, pairs[:, 1], pairs[:, 0])
        second = np.where(swap, pairs[:, 0], pairs[:, 1])
        null[index] = float(np.mean(np.where(first > second, 1.0, np.where(first == second, 0.5, 0.0))))
    return {
        "observed_pair_accuracy": observed,
        "shuffle_mean": float(null.mean()),
        "shuffle_ci_low": float(np.quantile(null, 0.025)),
        "shuffle_ci_high": float(np.quantile(null, 0.975)),
        "one_sided_p_observed_or_higher": float((1 + np.sum(null >= observed)) / (1 + len(null))),
        "n_permutations": int(n_permutations),
        "resampling_unit": "mutation; independently swap the two partner scores",
    }
