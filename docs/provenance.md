# Calculation provenance

This release contains a small reproduction subset of EdgeLossAudit, not the full upstream training or report-generation package.

- `base_repo/directional_esignet/intact_mp_robustness.py`: per-stratum statistics and mutation-cluster inference.
- `base_repo/directional_esignet/intact_mp_multiplicity.py`: paired inference and multiplicity correction.
- `base_repo/directional_esignet/ppi_rewirebench/brca1_brct_dscript_audit.py`: BRCA1-BRCT partner-score shuffle controls.

`statistics.py` retains the original function bodies for per-stratum statistics, paired mutation-cluster inference, BH/Holm correction and BRCA1-BRCT shuffling. Training, sequence retrieval and feature generation are outside this release.
