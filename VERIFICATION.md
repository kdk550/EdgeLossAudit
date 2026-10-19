# Recorded core-evidence verification

Date: 2026-10-19. Data release: `core-evidence-20261019`.

The source tree and its dataset were copied outside the project workspace, installed into a separate Python target directory with `pip install --no-deps --no-build-isolation`, and run from that independent directory. That run used the full local dataset; the public repository carries the code and dataset documentation only. Runtime dependencies were the already installed versions recorded in `requirements-tested.txt`. This checks packaging and independence from the study workspace; it does not claim a new operating system or dependency solver was tested.

Command after installation:

```bash
python -m edge_loss_audit --data dataset --out results --full-inference
```

Result: **passed**. The combined data set, containing 31 dataset files, passed a fresh independent installation and complete run in **75.30 seconds** on the current host CPU. This is an observed local runtime, not a performance guarantee. No separate dataset download is required. The command checks scientific reference values, cohort definitions, seed coverage and fold consistency. The data files themselves are not redistributed in this repository; see `README.md` for their sources.

| Check | Result |
| --- | --- |
| Primary common cohort and all 17 comparator AP/AUROC/conditional metrics | Matches retained references |
| Primary WT-existing target before common-row restriction | 1,613 rows, 586 positives |
| Primary masked marginal | AP 0.522257; AUROC 0.645488; conditional macro AUROC and pair accuracy 0.500000 |
| Component pair-aware conditional metrics | Macro AUROC 0.592770; pair accuracy 0.625966 |
| Homology-proxy pair-mean probe | Conditional macro AUROC 0.541468 |
| UniRef50 identity-closed delta-site probe | Pooled AP 0.425222; conditional macro AUROC 0.507174; a different model/split from the preceding row |
| Main IntAct seed-averaged stratum statistics | Reconstructed from all ten released row-score replicates |
| Complete IntAct inference family | All 48 point estimates, bootstrap intervals, group/pair counts and sign permutations verified |
| IntAct publication-matched primary pair-accuracy difference | -0.0307336; 95% CI [-0.0623007, 0.0067423] |
| IntAct secondary macro-AUROC correction | BH 0.0079996; Holm 0.0599970 |
| IGVF conditional pair-accuracy permutation probability | 0.5841584, from 100 saved permutation metrics |
| BRCA1-BRCT complete-loss / sensitivity pair-accuracy shuffle probability | 0.1764912 / 0.0952952; both regenerated with 20,000 swaps |
| Known-complex candidate coverage | 84/1,533 observations; 24/553 positives |

The inference comparison permits at most two Monte Carlo counts of floating-point tie sensitivity in a 20,000-draw sign-permutation probability. Point estimates and bootstrap intervals are compared with absolute tolerance 1e-10. BH/Holm corrections of the retained 48 reference P values are also independently recalculated.

The actual run generated `metrics.csv`, `intact_inference.csv`, `multiplicity_check.csv` and `summary.json` in a temporary results directory; generated outputs are not duplicated in this compact repository. No model fitting, external downloads, primary or BRCA1-BRCT confidence-interval regeneration, full sequence identity recomputation or figure rendering was performed. Non-main IntAct stratum statistics are sufficient to reproduce the provided inference, but their omitted observation-level predictions cannot be independently reconstructed from this release.
