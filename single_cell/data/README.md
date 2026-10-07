# Input data

Raw sequencing data for this neutrophil single-cell RNA-seq experiment are deposited in the DDBJ Sequence Read Archive (DRA) under accession **DRA031926**.

Processed per-sample Cell Ranger filtered feature-barcode matrices are deposited in the DDBJ Genomic Expression Archive (GEA) under accession **E-GEAD-1320**.

For the public Seurat workflow, place the following four downloaded HDF5 files directly under `Rawdata/`:

```text
Rawdata/
  PBS_1_sample_filtered_feature_bc_matrix.h5
  PBS_2_sample_filtered_feature_bc_matrix.h5
  Alb_1_sample_filtered_feature_bc_matrix.h5
  Alb_2_sample_filtered_feature_bc_matrix.h5
```

The public scripts use sample labels `PBS-1`, `PBS-2`, `Alb-1`, and `Alb-2` internally while reading the deposited files shown above.

Large intermediate Seurat objects (`*.rds`) are intentionally not distributed in this repository. They are generated sequentially by the scripts under `scripts/scRNA/` and written to `rds/`.
