# EdgeLossBench / EdgeLossAudit

Code for checking whether missense mutations preferentially break particular protein-protein interactions. Data and predictor code are not in this repo; see [Where the data comes from](#where-the-data-comes-from).

## Background

Two protein may keep interacting after a mutation, stop interacting, or interact in a new way. A model that only reports one pooled score over all mutation-partner pairs hides which of these it can actually tell apart, so we went back over several published PPI resources and recomputed the stratified metrics, the identity controls and the multiple-testing corrections ourselves.

The repository runs the same checks that we used while writing the analysis. Every number it prints is compared to a value recorded in the study's reference tables. If anything drifts, the script stops with `Reference mismatch` instead of printing a result you might trust.

## Running it

Python 3.10 or newer; the last run was on 3.12.3. Install into a fresh environment:

```bash
pip install -r requirements-tested.txt
pip install --no-deps .
```

Then point `--data` at a directory holding the five cohort folders (`primary/`, `intact_mp/`, `igvf/`, `brca1_brct/`, `structure/`):

```bash
python -m edge_loss_audit --data /path/to/data --out results --full-inference
```

Without `--full-inference` the script checks the point estimates and the multiplicity corrections, which takes a few seconds. With it, it also runs the 5,000 mutation-cluster bootstrap draws and 20,000 sign permutations for all 48 IntAct comparisons, and the BRCA1-BRCT shuffle tests. That was 75 s on the machine we used.

Outputs are `metrics.csv`, `intact_inference.csv`, `multiplicity_check.csv` and `summary.json` in `--out`. What was checked on the last run, and what was not, is in [VERIFICATION.md](VERIFICATION.md).

## What is checked

| Claim | What the script recomputes |
| --- | --- |
| A pooled ranking does not show that a model discriminates between partners | AP and AUROC over all 1,533 common observations, plus the mutation-conditional AUROC and pair accuracy, for all 17 primary comparators |
| The conditional margin shrinks once identity is controlled for | The primary identity-closed scores, and the IntAct-MP identity30 folds and ten seed replicates |
| The IntAct comparison matched to the original publication does not show an improvement | Seed-averaged per-stratum statistics, mutation-cluster inference and BH/Holm over all 48 tests |
| The secondary signal depends on the endpoint and the correction used | All six label scenarios are kept in the family; no favourable tests are dropped |
| The independent assays stay separate | IGVF pair-accuracy permutation probability from 100 saved permutation runs; BRCA1-BRCT D-SCRIPT shuffle controls at 20,000 swaps |
| Structural support is thin | Known-complex coverage among the 1,533 observations |

## Where the data comes from

The cohort files are converted versions of public resources. They are not in this repository, mostly because they are derived from someone else's data and we would rather people pull the originals from the source.

| Folder | Source | Link |
| --- | --- | --- |
| `primary/` | eSIG-Net (Zenodo, CC BY 4.0). We retrained the reference discriminator on our labels and splits. | <https://zenodo.org/records/15486737> |
| `intact_mp/` | IntAct mutation-impact annotations (CC BY 4.0). Interaction accessions and PMIDs are kept. | <https://www.ebi.ac.uk/intact/> |
| `igvf/` | IGVF semi-qY2H edge-perturbation scores. We used file `IGVFFI1611ZIBI` from analysis set `IGVFDS9876GIEQ`. | <https://data.igvf.org/tabular-files/IGVFFI1611ZIBI/> |
| `brca1_brct/` | Mangkusaputra et al. 2025, *Nucleic Acids Res*, doi:10.1093/nar/gkaf848. Dataset on Mendeley Data (CC BY 4.0). | <https://data.mendeley.com/datasets/9fckzvn6wf/1> |
| `structure/` | UniProt accessions matched against shared PDB entries. The match was done by us from the UniProt and PDB releases current at the time. | <https://www.uniprot.org>, <https://www.rcsb.org> |

We did not record the IntAct download date, which we should go back and fix before this is cited anywhere. The IGVF file `IGVFFI1611ZIBI` was released on 2026-02-13 and has since been replaced by `IGVFFI2460BBXY` in the IGVF portal; all IGVF numbers here come from the February file.

The file layout, column meanings and the processing applied to each upstream file are in [dataset/README.md](dataset/README.md), `dataset/DATA_DICTIONARY.tsv`, `dataset/PROVENANCE.tsv` and `dataset/MANIFEST.tsv`.

## Predictor code

None of the predictor code is vendored here. Frozen D-SCRIPT / Topsy-Turvy scores use the `samsl/topsy_turvy_human_v1` weights from the authors' repository (<https://github.com/samsledje/D-SCRIPT>). eSIG-Net is from <https://zenodo.org/records/15486737>. IGVF and IntAct are accessed through their own portals, linked above.

## Not in this release

This is the subset we could package up for a reviewer to rerun; it is not the full study pipeline.

- Databases are not downloaded and embeddings are not extracted.
- No models are trained or refit, apart from what was already done upstream. The IGVF probes are not rerun; the permutation probability is recomputed from the 100 saved permutation metrics.
- Confidence intervals for the primary benchmark and for BRCA1-BRCT are kept as reported, not regenerated.
- IntAct row predictions for the non-main label scenarios are omitted. Their seed-averaged per-stratum statistics are included, which is enough to rebuild the 48-test family.
- Figure scripts are elsewhere.

## Licence

The code in this repository is MIT (see [LICENSE](LICENSE)). The data keeps whatever licence its source has; the MIT licence does not cover it.
