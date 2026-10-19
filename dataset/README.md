# Core evidence for EdgeLossBench / EdgeLossAudit

Release identifier: `core-evidence-20261019`. This file describes the data layout expected by `edge_loss_audit`. The data files are **not redistributed**: only this documentation and the data dictionary are kept in the repository, and the cohort directories are obtained from the upstream sources below. This is a selected evidence release, distinct from the earlier six-track PPI-RewireBench v0.2 dataset. It includes no manuscript, rendered figures, model weights or embedding caches.

## Contents and endpoints

| Directory | Observations and endpoint | Role |
| --- | --- | --- |
| `primary/` | 1,533 common observations, 553 edge-loss positives, 462 mutations and 17 comparator scores; the input germline observations and identity-closed predictions are also included | Primary WT-existing binary edge-loss benchmark |
| `intact_mp/` | 36,158 curated observations (12,112 positives); 35,861 sequence-valid rows; 7,138 main-label mixed-outcome rows for each of ten identity30 seed replicates | Independent curated interaction-impact test |
| `igvf/` | 1,839 released assay observations; 1,833 affected-sequence-valid mutation-score rows; 1,827 rows with partner-proxy scores | Quantitative semi-qY2H, binary reporting view `log2fc <= -1` |
| `brca1_brct/` | 8,128 pair-complete observations from 4,064 mutations; 542 frozen D-SCRIPT rows from the union of the two mixed-outcome cohorts | Complete-loss threshold `<= -5`; separate sensitivity threshold `<= -1` |
| `structure/` | Known-complex support flags for the 1,533 common observations | Coverage only: 84 covered rows, 24 covered positives |

### Primary benchmark

`observations.csv.gz` preserves all columns of the 17-comparator common-row score table. Two distinct identity-control protocols are supplied:

- `homology_proxy_predictions.csv.gz`: five held-out folds of the homology-proxy logistic pair-mean probe; this is the source of conditional macro AUROC 0.5415.
- `uniref50_identity_predictions.csv.gz`: five held-out folds of the UniRef50-seeded identity-closed logistic delta-site probe; this is the source of pooled AP 0.4252. Its conditional macro AUROC is 0.5072, a different protocol/model from the 0.5415 result.

These identity controls do not establish exclusion from language-model pretraining. `uniref50_identity_audit.json` is the original audit summary, not a new sequence-similarity calculation. The complete sequence files and train/validation construction workflow are not included here.

`raw_observations.csv.gz` names the columns of the original headerless germline file explicitly: the first identifier is the mutant, `mutant_label` is its interaction state, and `wt_label` is the wild-type interaction state. Restrict to `wt_label == 1` and set loss to `mutant_label == 0`. This yields 1,613 observations and 586 positives before common-row restriction. Reversing the biological universe is different from reversing a model score.

`reported_comparators.csv`, `direction_reference.json` and `statistical_reference.json` preserve reported comparison and uncertainty references. The reproduction command recalculates point metrics; it does not regenerate the primary benchmark confidence intervals.

### IntAct-MP and multiplicity

Positives are decreasing/disrupting annotations; negatives are no-effect/increasing annotations. These are curated impact labels rather than direct WT-existing binary edge-loss measurements. The identity30 split controls the affected-protein side, with minimum alignment coverage 0.8. Partner-side cluster membership is provided for audit and is not a promise of dual-side isolation.

`main_mixed_seed_predictions.csv.gz` contains all mixed-outcome observation scores for all ten main-label identity30 fits. Mixed outcomes are selected by labels, independently of model performance. Pure-outcome mutations are excluded from this conditional-score file and remain in `observations.csv.gz`.

`group_statistics.csv.gz` contains per-stratum conditional statistics averaged across the ten seeds, for all six label scenarios, two grouping rules and three models. Averaging is performed on per-seed stratum statistics, not on raw scores. Its fields allow reconstruction of the complete 48-test inference family. The reproduction command also verifies the main-label sufficient statistics against the ten released row-score replicates. Other scenarios' row predictions are omitted.

Publication matching compares positive and negative partners within mutation-by-PubMed strata, then resamples mutation clusters. Macro AUROC weights mutation clusters equally; micro pair accuracy weights positive-negative pairs. `inference_protocol.json` fixes the primary endpoint, random seed and resampling counts. `reported_inference.csv` retains all 48 original tests and correction results; `seed_metrics.csv` retains the broader per-seed summaries, including split comparisons whose full-cohort row scores are not in this compact release.

### Independent tracks

IGVF partner sequences are canonical translation proxies, not verified assay ORF sequences. The reproduction command recalculates the one-sided permutation probability from the 100 saved permutation metrics; it does not refit those probes or reconstruct the shuffled datasets.

BRCA1-BRCT is restricted to two partners of one domain. Two mutations with exact IntAct mutation-partner overlap were removed as complete pairs. Its D-SCRIPT score is `P(WT, partner) - P(mutant, partner)` from the frozen pretrained model `samsl/topsy_turvy_human_v1`. D-SCRIPT coverage is complete for the 542-row conditional union, not for the full 8,128-row cohort. The two endpoint-specific shuffle tests use 20,000 swaps with seeds 20360720 and 20361720 respectively. `cohort_protocol.json` records cohort filtering; `pair_reference.json` retains the original endpoint statistics.

## Sources, provenance and reuse

- Primary benchmark: the study's processed derivative of the eSIG-Net resource (CC BY 4.0), <https://zenodo.org/records/15486737>. The reference discriminator was retrained for this study's labels and splits.
- IntAct-MP: processed IntAct mutation-impact annotations (CC BY 4.0), <https://www.ebi.ac.uk/intact/>. Interaction accessions and publication identifiers are retained. The exact database release/download date is not established in this compact inventory and should be supplied before final deposition.
- IGVF: file `IGVFFI1611ZIBI`, dataset `IGVFDS9876GIEQ`, <https://data.igvf.org/tabular-files/IGVFFI1611ZIBI/>. This file was released 2026-02-13 and later superseded by `IGVFFI2460BBXY` (`1percent_semiqY2H_edge_scores.tsv.gz`); the analysis used the 2026-02-13 version.
- BRCA1-BRCT: Mangkusaputra et al. (2025), article DOI `10.1093/nar/gkaf848`; original dataset DOI `10.17632/9fckzvn6wf.1` (CC BY 4.0), <https://data.mendeley.com/datasets/9fckzvn6wf/1>.
- Structure flags are derived from the study's UniProt/shared-PDB mapping (<https://www.uniprot.org/>, <https://www.rcsb.org/>). Underlying sequence and structure download caches are omitted.
- Method code: D-SCRIPT / Topsy-Turvy, <https://github.com/samsledje/D-SCRIPT> (MIT); frozen weights `samsl/topsy_turvy_human_v1` are not redistributed.

`PROVENANCE.tsv` records each asset's project-relative source path, destination and transformation. The paths describe the study workspace; they are not download URLs or runtime dependencies. `MANIFEST.tsv` lists released file sizes and row counts. `DATA_DICTIONARY.tsv` describes CSV/TSV fields. Missing numeric values are empty CSV fields; some source metadata retain the literal `nan` or an empty string. Model scores are generally unitless and are not assumed calibrated probabilities.

No new license is assigned to these data, and no data files are redistributed in this repository. Anyone reproducing them must follow the upstream resource terms and attribution requirements listed above; author and archive citation metadata are not yet supplied.
