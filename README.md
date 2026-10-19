# EdgeLossBench / EdgeLossAudit

Code for checking whether missense mutations preferentially break particular protein-protein interactions. Data and predictor code are not in this repo; see [Input data sources](#input-data-sources).

## Motivation

Two protein may keep interacting after a mutation, stop interacting, or interact in a new way. A model that only reports one pooled score over all mutation-partner pairs hides which of these it can actually tell apart, so we went back over several published PPI resources and recomputed the stratified metrics, the identity controls and the multiple-testing corrections ourselves.

The repository runs the same checks that we used while writing the analysis. Every number it prints is compared to a value recorded in the study's reference tables. If anything drifts, the script stops with `Reference mismatch` instead of printing a result you might trust.

## Usage

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

## What the script checks

| Claim | What the script recomputes |
| --- | --- |
| Partner discrimination beyond pooled ranking | AP and AUROC over all 1,533 common observations, plus the mutation-conditional AUROC and pair accuracy, for all 17 primary comparators |
| Effect of identity control | The primary identity-closed scores, and the IntAct-MP identity30 folds and ten seed replicates |
| IntAct publication-matched comparison | Seed-averaged per-stratum statistics, mutation-cluster inference and BH/Holm over all 48 tests |
| Endpoint and correction dependence | All six label scenarios are kept in the family |
| Independent assays | IGVF pair-accuracy permutation probability from 100 saved permutation runs; BRCA1-BRCT D-SCRIPT shuffle controls at 20,000 swaps |
| Structural coverage | Known-complex coverage among the 1,533 observations |

## Input data sources

The cohort files are derived from the public resources below, and the originals should be taken from the sources themselves.

| Folder | Source | Link |
| --- | --- | --- |
| `primary/` | eSIG-Net (Zenodo, CC BY 4.0). The reference discriminator was retrained on our labels and splits. | <https://zenodo.org/records/15486737> |
| `intact_mp/` | IntAct mutation-impact annotations (CC BY 4.0), with interaction accessions and PMIDs kept. | <https://www.ebi.ac.uk/intact/> |
| `igvf/` | IGVF semi-qY2H edge-perturbation scores, file `IGVFFI1611ZIBI` (released 2026-02-13) from analysis set `IGVFDS9876GIEQ`. | <https://data.igvf.org/tabular-files/IGVFFI1611ZIBI/> |
| `brca1_brct/` | Mangkusaputra et al. 2025, *Nucleic Acids Res*, doi:10.1093/nar/gkaf848. Dataset on Mendeley Data (CC BY 4.0). | <https://data.mendeley.com/datasets/9fckzvn6wf/1> |
| `structure/` | Protein accessions mapped to UniProt (<https://www.uniprot.org>) and shared entries in PDB (<https://www.rcsb.org>); the mapping is computed in the study from these releases. | <https://www.uniprot.org>, <https://www.rcsb.org> |

The file layout, column meanings and the processing applied to each upstream file are in [dataset/README.md](dataset/README.md), `dataset/DATA_DICTIONARY.tsv`, `dataset/PROVENANCE.tsv` and `dataset/MANIFEST.tsv`.

## Predictor sources

Predictors are used from the authors' own releases: D-SCRIPT / Topsy-Turvy with the frozen `samsl/topsy_turvy_human_v1` weights (<https://github.com/samsledje/D-SCRIPT>) and eSIG-Net (<https://zenodo.org/records/15486737>). IntAct and IGVF are read through their portals linked above.

## Licence

The code in this repository is MIT (see [LICENSE](LICENSE)). The data keeps whatever licence its source has; the MIT licence does not cover it.
