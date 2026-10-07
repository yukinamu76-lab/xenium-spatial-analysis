# Public-release audit

## Checks completed

- No local absolute paths such as `/Users/...` or `/Volumes/...` were found in the supplied analysis-code tree.
- No email addresses, submitter identifiers, GEA submission IDs, BioSample IDs, or other internal submission metadata were found in the supplied analysis-code tree.
- Public input paths in `scripts/scRNA/00_setup.R` match the four per-sample GEA HDF5 filenames.
- Raw-data accession **DRA031926** and processed-data GEA accession **E-GEAD-1320** are documented in the root README and `data/README.md`.
- Separate reproducibility records are included for the two analysis environments:
  - `sessionInfo_scRNA.txt` for the scRNA-seq workflow (R 4.6.0)
  - `sessionInfo_GSEA.txt` for the GSEA workflow (R 4.6.1)
- Environment-capture scripts are retained in `scripts/reproducibility/`.

## Recommended final verification

Before creating the archival release, rerun the workflow from the deposited HDF5 matrices in the corresponding analysis environments and confirm that the reported cluster assignments, Fig. 6/EV4 values, and exported tables agree with the submitted manuscript.
