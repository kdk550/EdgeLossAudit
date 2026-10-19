"""Recalculate the paper's core evidence from the compact data release."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from threadpoolctl import threadpool_limits

from .statistics import adjust_bh, adjust_holm, paired_inference, per_group_statistics, shuffle_pair_accuracy


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False, float_precision="round_trip")


def check_close(actual, expected, label, atol=1e-10):
    if not np.allclose(actual, expected, rtol=0, atol=atol, equal_nan=True):
        raise ValueError(f"Reference mismatch: {label}")


def group_stats(frame, score, columns=("mutation",), target="target"):
    work = frame.rename(columns={target: "analysis_target"}) if target != "analysis_target" else frame
    mixed = work.groupby(list(columns), sort=False, dropna=False)["analysis_target"].transform("nunique").eq(2)
    return per_group_statistics(work.loc[mixed], score, list(columns))


def metrics(frame, score, cohort, target="target"):
    y = frame[target].to_numpy(dtype=int)
    values = frame[score].to_numpy(dtype=float)
    if not set(np.unique(y)).issubset({0, 1}) or not np.isfinite(values).all():
        raise ValueError(f"Invalid labels or scores: {cohort}/{score}")
    groups = group_stats(frame, score, target=target)
    return {
        "cohort": cohort, "method": score, "rows": len(frame), "positives": int(y.sum()),
        "prevalence": float(y.mean()), "average_precision": float(average_precision_score(y, values)),
        "auroc": float(roc_auc_score(y, values)), "mixed_mutations": len(groups),
        "within_macro_auroc": float(groups.macro_auroc.mean()),
        "within_pair_accuracy": float(groups.pair_correct.sum() / groups.pairs.sum()),
        "positive_negative_pairs": int(groups.pairs.sum()),
        "score_constant_within_mutation": bool(frame.groupby("mutation")[score].nunique().eq(1).all()),
    }


def validate_main_seed_statistics(predictions, released_stats):
    """Check the compact sufficient statistics against all ten released row-score replicates."""
    methods = {"pair_aware": "pair_aware_linear", "mutation_only": "mutation_only_linear", "partner_shuffle": "partner_shuffled_linear"}
    for grouping, columns in [("mutation", ["mutation"]), ("mutation_pubmed", ["mutation", "pubmed_id"])]:
        for method, score in methods.items():
            frames = []
            for _, frame in predictions.groupby("seed", sort=False):
                frames.append(group_stats(frame, score, columns, "analysis_target"))
            averaged = pd.concat(frames).groupby("group", as_index=False).agg(
                macro_auroc=("macro_auroc", "mean"), pair_correct=("pair_correct", "mean"), pairs=("pairs", "first"))
            expected = released_stats.query("scenario == 'main_all' and grouping == @grouping and method == @method")
            merged = averaged.merge(expected, on="group", suffixes=("_new", "_released"), validate="one_to_one")
            if len(merged) != len(averaged) or len(merged) != len(expected):
                raise ValueError("Missing main-label group statistics")
            for field in ["macro_auroc", "pair_correct", "pairs"]:
                check_close(merged[field + "_new"], merged[field + "_released"], f"main seed statistics/{grouping}/{method}/{field}")


def inference_from_groups(stats, protocol, resample):
    rows = []
    comparisons = {"pair_aware_vs_mutation_only": ("pair_aware", "mutation_only"),
                   "pair_aware_vs_partner_shuffle": ("pair_aware", "partner_shuffle")}
    for si, scenario in enumerate(sorted(stats.scenario.unique())):
        for gi, grouping in enumerate(["mutation", "mutation_pubmed"]):
            selected = stats.query("scenario == @scenario and grouping == @grouping")
            for ci, (comparison, (first, second)) in enumerate(comparisons.items()):
                left = selected[selected.method.eq(first)]
                right = selected[selected.method.eq(second)]
                for mi, metric in enumerate(["macro_auroc", "pair_accuracy"]):
                    if resample:
                        result = paired_inference(left, right, metric, protocol["n_bootstrap"], protocol["n_permutations"],
                                                  protocol["seed"] + si * 1000 + gi * 100 + ci * 10 + mi)
                    else:
                        merged = left.merge(right, on="group", suffixes=("_first", "_second"), validate="one_to_one")
                        if metric == "macro_auroc":
                            delta = merged.macro_auroc_first - merged.macro_auroc_second
                            result = {"delta": float(delta.groupby(merged.resampling_cluster_first, sort=False).mean().mean())}
                        else:
                            result = {"delta": float((merged.pair_correct_first - merged.pair_correct_second).sum() / merged.pairs_first.sum())}
                    rows.append({"scenario": scenario, "split": protocol["split"], "grouping": grouping,
                                 "comparison": comparison, "metric": metric, **result})
        if resample:
            print(f"Recomputed inference: {scenario}", flush=True)
    return pd.DataFrame(rows)


def run(data, out, full_inference):
    out.mkdir(parents=True, exist_ok=True)
    results = []
    primary = read(data / "primary/observations.csv.gz")
    reference = read(data / "primary/reported_comparators.csv")
    for row in reference.itertuples(index=False):
        score_column = "prevalence_null_score" if row.method == "prevalence_null" else row.method
        result = metrics(primary, score_column, "primary_component_common")
        for field in ["average_precision", "auroc"]:
            check_close(result[field], getattr(row, field), f"primary/{row.method}/{field}")
        check_close(result["within_macro_auroc"], row.within_macro_auroc, f"primary/{row.method}/macro")
        check_close(result["within_pair_accuracy"], row.within_micro_pair_accuracy, f"primary/{row.method}/pairs")
        results.append(result)
    if len(primary) != 1533 or primary.target.sum() != 553:
        raise ValueError("Primary common cohort changed")
    raw = read(data / "primary/raw_observations.csv.gz")
    wt = raw[raw.wt_label.eq(1)]
    if len(wt) != 1613 or wt.mutant_label.eq(0).sum() != 586:
        raise ValueError("Primary biological direction changed")
    homology = read(data / "primary/homology_proxy_predictions.csv.gz")
    for method in ["masked_marginal", "pair_aware_probe"]:
        result = metrics(homology, method, "primary_homology_proxy")
        if method == "pair_aware_probe":
            check_close(result["within_macro_auroc"], 0.5414678991318336, "homology proxy conditional AUROC")
        results.append(result)
    identity = read(data / "primary/uniref50_identity_predictions.csv.gz")
    for method in ["masked_marginal", "pair_aware_probe"]:
        result = metrics(identity, method, "primary_uniref50_identity")
        if method == "pair_aware_probe":
            check_close(result["average_precision"], 0.425221933502125, "UniRef50 identity pooled AP")
        results.append(result)

    stats = read(data / "intact_mp/group_statistics.csv.gz")
    seed_rows = read(data / "intact_mp/main_mixed_seed_predictions.csv.gz")
    validate_main_seed_statistics(seed_rows, stats)
    protocol = json.loads((data / "intact_mp/inference_protocol.json").read_text())
    if sorted(seed_rows.seed.unique()) != sorted(protocol["seeds"]):
        raise ValueError("Incomplete seed set")
    cohort = read(data / "intact_mp/observations.csv.gz")
    if len(cohort) != 36158 or cohort.target.sum() != 12112 or cohort.sequence_valid.sum() != 35861:
        raise ValueError("IntAct cohort changed")
    fold_counts = cohort.groupby("affected_identity30_cluster").identity30_fold.nunique()
    if fold_counts.max() != 1:
        raise ValueError("Affected-protein identity clusters cross folds")
    # Full correction family is mandatory even though only main labels are emphasized.
    reported = read(data / "intact_mp/reported_inference.csv")
    if len(reported) != 48 or reported.primary_endpoint.sum() != 1:
        raise ValueError("The fixed 48-test correction family changed")
    reported["recomputed_bh_all"] = adjust_bh(reported.permutation_p.to_numpy())
    reported["recomputed_holm_all"] = adjust_holm(reported.permutation_p.to_numpy())
    secondary = ~reported.primary_endpoint
    check_close(reported.recomputed_bh_all, reported.p_bh_all, "BH correction")
    check_close(reported.recomputed_holm_all, reported.p_holm_all, "Holm correction")
    check_close(adjust_bh(reported.loc[secondary, "permutation_p"].to_numpy()), reported.loc[secondary, "p_bh_secondary"], "secondary BH")
    inference = inference_from_groups(stats, protocol, full_inference)
    keys = ["scenario", "split", "grouping", "comparison", "metric"]
    joined = inference.merge(reported, on=keys, suffixes=("_new", "_reported"), validate="one_to_one")
    check_close(joined.delta_new, joined.delta_reported, "all 48 inference point estimates")
    if full_inference:
        for field in ["ci_low", "ci_high", "groups", "matched_strata", "positive_negative_pairs", "probability_gt_zero"]:
            check_close(joined[field + "_new"], joined[field + "_reported"], f"resampling/{field}")
        # A null statistic can tie the observation at floating-point precision.
        check_close(joined.permutation_p_new, joined.permutation_p_reported, "sign permutations", atol=2 / (protocol["n_permutations"] + 1))
        inference["p_bh_all"] = adjust_bh(inference.permutation_p.to_numpy())
        inference["p_holm_all"] = adjust_holm(inference.permutation_p.to_numpy())
    inference.to_csv(out / "intact_inference.csv", index=False)
    reported.to_csv(out / "multiplicity_check.csv", index=False)

    mutation_igvf = read(data / "igvf/mutation_scores.csv.gz")
    results.append(metrics(mutation_igvf, "masked_marginal_loss_score", "igvf_affected_sequence_valid"))
    pair_igvf = read(data / "igvf/pair_scores.csv.gz")
    igvf = metrics(pair_igvf, "pair_aware_oof_score", "igvf_partner_proxy")
    results.append(igvf)
    shuffle = read(data / "igvf/shuffle_metrics.csv")
    igvf_p = float((1 + shuffle.within_mutation_pair_accuracy.ge(igvf["within_pair_accuracy"]).sum()) / (1 + len(shuffle)))
    expected = json.loads((data / "igvf/shuffle_reference.json").read_text())
    check_close(igvf_p, expected["one_sided_permutation_p"]["within_mutation_pair_accuracy_observed_greater_than_shuffle"], "IGVF shuffle P")

    mutation_brct = read(data / "brca1_brct/mutation_scores.csv.gz")
    pair_brct = read(data / "brca1_brct/pair_scores.csv.gz")
    brct_reference = json.loads((data / "brca1_brct/pair_reference.json").read_text())
    brct_controls = {}
    for target, endpoint, seed in [
        ("target_complete_loss", "source_defined_complete_loss_score_le_neg5", 20360720),
        ("target_neg1_sensitivity", "sensitivity_score_le_neg1", 20361720),
    ]:
        results.append(metrics(mutation_brct, "masked_marginal_loss_score", f"brca1_brct/{target}", target))
        result = metrics(pair_brct, "dscript_delta_loss_score", f"brca1_brct_dscript/{target}/conditional_union", target)
        # The supplied D-SCRIPT cohort is a union of mixed mutations from the two endpoints.
        # Its pooled metrics do not describe the full 8128-row assay cohort.
        results.append(result)
        control = shuffle_pair_accuracy(pair_brct.assign(target=pair_brct[target]), "dscript_delta_loss_score", 20000, seed)
        expected_control = brct_reference["endpoint_results"][endpoint]["partner_shuffle_control"]
        check_close(control["one_sided_p_observed_or_higher"], expected_control["one_sided_p_observed_or_higher"], f"BRCT shuffle/{target}")
        brct_controls[target] = control

    coverage = read(data / "structure/coverage.csv.gz")
    covered = coverage.known_complex_candidate.astype(bool)
    support = {"rows": len(coverage), "known_complex_rows": int(covered.sum()),
               "positives": int(coverage.lost_target.sum()), "known_complex_positives": int(coverage.loc[covered, "lost_target"].sum())}
    if (support["known_complex_rows"], support["known_complex_positives"]) != (84, 24):
        raise ValueError("Known-complex coverage changed")
    pd.DataFrame(results).to_csv(out / "metrics.csv", index=False)
    primary_endpoint = reported.loc[reported.primary_endpoint].iloc[0]
    secondary_endpoint = reported.query("scenario == 'main_all' and grouping == 'mutation' and comparison == 'pair_aware_vs_partner_shuffle' and metric == 'macro_auroc'").iloc[0]
    summary = {
        "status": "passed", "full_inference_resampling": full_inference,
        "primary_common_rows": len(primary), "primary_methods_recomputed": len(reference),
        "intact_seed_count": len(protocol["seeds"]), "intact_primary_delta": float(primary_endpoint.delta),
        "intact_primary_ci": [float(primary_endpoint.ci_low), float(primary_endpoint.ci_high)],
        "intact_secondary_bh": float(secondary_endpoint.p_bh_all), "intact_secondary_holm": float(secondary_endpoint.p_holm_all),
        "igvf_one_sided_shuffle_p": igvf_p, "brca1_brct_controls": brct_controls, "structure_support": support,
        "limits": ["No upstream sequence retrieval, model fitting, scoring or figure generation is performed.",
                   "Non-main IntAct scenarios are supplied as seed-averaged group statistics; their row predictions are omitted.",
                   "IGVF permutation P is recalculated from the 100 saved permutation metrics, without refitting probes.",
                   "Primary benchmark and BRCA1-BRCT confidence intervals are not regenerated by this command."]}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--full-inference", action="store_true", help="Rerun all 48 IntAct cluster bootstrap/sign-permutation tests; otherwise verify points and multiplicity corrections.")
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        run(args.data.resolve(), args.out.resolve(), args.full_inference)


if __name__ == "__main__":
    main()
