# Code–Methods audit

## Confirmed alignment

- QC: `min.cells = 3`; `200 < nFeature_RNA < 6000`; `percent.mt < 15%`; mitochondrial genes matched with `^mt-`.
- All-cell and neutrophil-only preprocessing: LogNormalize, 2,000 variable genes (`vst`), regression of `percent.mt`, PCA with 30 components, Harmony using `sample_id` to obtain a shared low-dimensional representation, Harmony dimensions 1–15, clustering resolution 0.6, and UMAP on the same Harmony dimensions.
- Initial annotation uses canonical immune-cell markers and retains the neutrophil clusters; monocyte, NK, pDC, and basophil clusters are excluded before neutrophil-specific reclustering.
- Downstream cluster-vs-rest rankings exclude clusters 9 and 10 both as target clusters and from the `rest` background. Only clusters 0–8 are used for comprehensive GSEA.
- Slingshot removes clusters 9 and 10 before trajectory inference, uses the 2D UMAP embedding, leaves start and end clusters unspecified, and exports cell-level and cluster-summary pseudotime values.
- Comprehensive GSEA uses `fgseaMultilevel` separately for Hallmark, GO BP, and Reactome collections; complete cluster 0–8 results correspond to Dataset EV7.
- Targeted Fig. 6C GSEA evaluates five pathways across clusters 2, 6, 7, and 8 and applies BH correction across all 20 raw P values.
- Fig. EV4E uses the supplied v3.3 pathway-selection and ordering pipeline. `functions/01_read_data.R` is part of the active execution path; the other helper scripts are retained as provenance/supporting code.

## Manuscript-item nomenclature used in this release

- Supplementary Fig. 6 → **Fig. EV4**
- Supplementary Table 12 → **Table EV6**
- Supplementary Table 13 → **Table EV7**
- Supplementary Table 14 → **Table EV8**
- Supplementary Table 15 → **Dataset EV7**

## Public-release normalization performed

- Local absolute paths and personal identifiers were removed.
- Sample labels were normalized to `PBS-1`, `PBS-2`, `Alb-1`, and `Alb-2`.
- The final first-round neutrophil-cluster assignment was encoded in `00_setup.R` for reproducibility.
- Internal source paths were updated to reflect the public repository directory structure.
- Output-directory capitalization/naming was harmonized so the workflow works on case-sensitive systems.
- `09b_fig6B_pooled_cluster_ratio.R` reproduces the pooled calculation used for Fig. 6B and is explicitly a public-release reproduction script rather than an original analysis-stage script.
- Manuscript-item names in comments/output labels were updated to Fig. EV4, Table EV6–EV8, and Dataset EV7 without changing the analysis logic or statistical parameters.

## Before depositing a frozen release

1. DDBJ accessions are recorded as **DRA031926** for raw sequencing data and **E-GEAD-1320** for processed matrices.
2. Run the workflow from a clean R environment using the deposited/downloaded matrices.
3. Separate scRNA-seq and GSEA session information is provided in `sessionInfo_scRNA.txt` and `sessionInfo_GSEA.txt`; the corresponding capture scripts are retained under `scripts/reproducibility/`.
4. Confirm that generated cluster labels and figure values match the submitted manuscript.
5. Release this directory under the MIT License of the parent repository.
