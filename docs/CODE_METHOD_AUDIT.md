# Code–Methods audit

## Confirmed alignment

- QC retains cells with `nFeature_Xenium > 3`.
- The analysis uses a 50,000-cell leverage-score sketch, 50 PCs computed, and the first 20 Harmony dimensions downstream.
- Harmony integration is grouped by sample (`orig.ident`), with clustering at resolution 0.5.
- Cell-type annotation uses canonical markers plus cluster-enriched genes; `FindAllMarkers` uses `only.pos = TRUE`, `min.pct = 0.25`, and `logfc.threshold = 0.25`.
- CD8+ T-cell z scores are computed across all annotated CD8+ T cells before group-wise averaging.
- PMN-MDSC-like and exhaustion-associated CD8+ T-cell definitions use nuclear-localized transcripts as described in Methods.
- Spatial regions use DBSCAN (`eps = 50 um`, `minPts = 100`) and exclude regions with fewer than 5,000 cells.
- Logistic regression uses distance per 100 um, with Wald confidence intervals/P values; permutation testing uses 1,000 label permutations.
- Table/figure references have been updated to Fig. EV6, Table EV9, Table EV10, and Table EV11.

## Public-release normalization

- The local absolute project path was removed. Set `XENIUM_PROJECT_ROOT`, or run from the repository root/`R` directory.
- Input Xenium directories are expected under `raw_data/`; outputs are written under `results/`.
- Analysis parameters and statistical logic were not changed during the renumbering/public-format update.

## Public-release status

- DDBJ GEA accession `E-GEAD-1315` (array design `A-GEAD-246`) has been added to the README.
- An MIT software license has been added.
- Local absolute paths are excluded from the public scripts.

## Checks before creating a frozen archival release

1. Run from a clean R environment using the deposited/downloaded Xenium outputs.
2. Retain the generated `sessionInfo.txt` files.
3. Confirm that figure/table values match the submitted manuscript.
