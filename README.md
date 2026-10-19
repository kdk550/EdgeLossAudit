# EdgeLossBench / EdgeLossAudit: core evidence

Reproduction code for the study of partner-specific protein interaction loss after missense mutations. This repository contains the analysis code only. Third-party datasets and the upstream method code they build on are not redistributed; download them from the sources below.

## Run

Use Python 3.10 or newer (recorded with Python 3.12.3).

```bash
python -m pip install -r requirements-tested.txt
python -m pip install --no-deps .
python -m edge_loss_audit --data <data-dir> --out results --full-inference
```

`<data-dir>` must contain the five cohort directories (`primary/`, `intact_mp/`, `igvf/`, `brca1_brct/`, `structure/`) described in [dataset/README.md](dataset/README.md). They are not bundled. `dataset/README.md`, `dataset/DATA_DICTIONARY.tsv`, `dataset/PROVENANCE.tsv` and `dataset/MANIFEST.tsv` describe the file layout, fields and upstream processing; the data themselves must be obtained from the sources below.

Omit `--full-inference` for a faster check of point estimates and corrections. Full mode additionally runs the 5,000 mutation-cluster bootstrap draws and 20,000 two-sided sign permutations for all 48 IntAct-MP comparisons. Both modes verify the ten main-label seed replicates and rerun the BRCA1-BRCT partner-score shuffle test.

## Data sources

| Cohort | Upstream resource | Reference | License |
| --- | --- | --- | --- |
| `primary/` | eSIG-Net (the study's processed derivative; reference discriminator retrained) | <https://zenodo.org/records/15486737> | CC BY 4.0 |
| `intact_mp/` | IntAct mutation-impact annotations | <https://www.ebi.ac.uk/intact/> | CC BY 4.0 |
| `igvf/` | IGVF semi-qY2H edge-perturbation scores, file `IGVFFI1611ZIBI`, analysis set `IGVFDS9876GIEQ` | <https://data.igvf.org/tabular-files/IGVFFI1611ZIBI/> | Open access (no data-use certificate) |
| `brca1_brct/` | Mangkusaputra et al. (2025), *Nucleic Acids Research*, DOI `10.1093/nar/gkaf848` | <https://data.mendeley.com/datasets/9fckzvn6wf/1> | CC BY 4.0 |
| `structure/` | UniProt / PDB complex mapping derived in the study | <https://www.uniprot.org/>, <https://www.rcsb.org/> | see source sites |

The IGVF file `IGVFFI1611ZIBI` (released 2026-02-13) was later superseded in the IGVF portal by `IGVFFI2460BBXY` (`1percent_semiqY2H_edge_scores.tsv.gz`, 2026-08-10). The analysis reported here used the 2026-02-13 version; its row values therefore correspond to that release.

Each cohort is the study's processed derivative of the upstream resource, not a mirror of it. Derivative files follow the upstream licenses above and must retain their attribution if reproduced.

## Method code

Predictors and tools are used through their published code and, where applicable, released model weights; none are vendored here.

| Method | Code |
| --- | --- |
| D-SCRIPT / Topsy-Turvy (frozen `samsl/topsy_turvy_human_v1` weights) | <https://github.com/samsledje/D-SCRIPT> |
| eSIG-Net | <https://zenodo.org/records/15486737> |
| IntAct data access | <https://www.ebi.ac.uk/intact/> |
| IGVF data access | <https://data.igvf.org/> |

## Evidence reproduced

| Claim | Data and calculation |
| --- | --- |
| Pooled ranking does not establish partner discrimination | All 1,533 common observations and all 17 primary comparator scores; recalculate pooled AP/AUROC and mutation-conditional metrics |
| Conditional evidence weakens under identity control | Primary identity-closed row scores; IntAct-MP identity30 labels, folds, cluster membership and ten main-label mixed-outcome seed replicates |
| The primary publication-matched IntAct comparison does not demonstrate an improvement | Recalculate seed-averaged stratum statistics, mutation-cluster inference and all 48-test BH/Holm corrections |
| Secondary signal depends on endpoint and correction | Keep the full correction family, including all six label scenarios; no selection of only favorable tests |
| Independent assays remain separate | IGVF quantitative/proxy results with 100 saved permutation metrics, and BRCA1-BRCT mutation scores plus frozen D-SCRIPT conditional controls |
| Structural interpretation has limited support | Known-complex coverage labels for all 1,533 common observations |

`metrics.csv`, `intact_inference.csv`, `multiplicity_check.csv` and `summary.json` are generated under the requested output directory. A failed data/reference check exits with an error. See `VERIFICATION.md` for the recorded run, and `docs/provenance.md` for calculation provenance.

## Interpretation and scope

Higher scores mean greater predicted loss or perturbation. A constant score across the partners of one mutation yields conditional accuracy 0.5 by construction. Pooled association does not establish a biological severity mechanism. Each assay has its own endpoint; do not pool their labels or combine their metrics into one leaderboard.

This is a reproduction subset of EdgeLossAudit. It does not download raw databases, extract embeddings, train predictors, rerun the IGVF probes or regenerate the manuscript figures. Primary benchmark confidence intervals are retained as reported references rather than regenerated; BRCA1-BRCT confidence intervals are also not regenerated. Non-main IntAct row predictions are omitted; their complete seed-averaged stratum statistics allow the reported inference to be recalculated. The full upstream software and model-building workflow are outside this release.

## Citation and rights

Author metadata, the publication citation, software license and public repository URL remain to be supplied. Third-party data keep their source attribution and applicable terms; see the Data sources table above. No new license is assigned to third-party data by this snapshot.
