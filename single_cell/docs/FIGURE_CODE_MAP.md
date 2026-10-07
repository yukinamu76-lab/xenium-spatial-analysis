# Figure/table-to-code map

| Manuscript item | Analysis / content | Main code |
|---|---|---|
| Fig. EV4A | All-cell UMAP and cell-type annotation | `scripts/scRNA/01_qc_load.R`–`05_neutrophil_annotation.R`; `scripts/scRNA/18_supp_all_cell_UMAPs.R` |
| Fig. 6A | Neutrophil reclustering UMAP | `scripts/scRNA/06_neutrophil_subset_pca_elbow.R`, `scripts/scRNA/07_neutrophil_clustering_umap.R` |
| Fig. EV4B | QC evidence for cluster 9 exclusion | `scripts/scRNA/cluster9_10_exclusion_suppfig.R` |
| Fig. EV4C | Lineage-marker evidence for cluster 10 exclusion | `scripts/scRNA/cluster9_10_exclusion_suppfig.R` |
| Fig. 6B | Alb/PBS pooled cluster-composition ratio | `scripts/scRNA/09b_fig6B_pooled_cluster_ratio.R` |
| Fig. EV4D | Per-sample composition of excluded clusters | `scripts/scRNA/09_neutrophil_cluster_composition.R` and manuscript plotting workflow |
| Fig. 6C | Targeted GSEA of five representative programs | `scripts/gsea/targeted/targeted_GSEA_clusters_2_6_7_8.R` |
| Fig. EV4E | Broader functional GSEA landscape | `scripts/gsea/heatmap/run_all_v3.R` and sourced pipeline files |
| Fig. 6D | Slingshot pseudotime / Lineage 1 | `scripts/scRNA/11_neutrophil_trajectory.R` |
| Table EV6 | Sample-level QC, neutrophil counts/proportions, and cluster composition | `scripts/scRNA/01_qc_load.R`, `05_neutrophil_annotation.R`, `09_neutrophil_cluster_composition.R` |
| Table EV7 | Marker panels used for cell-type annotation and subcluster-purity assessment | marker definitions in `scripts/scRNA/00_setup.R`, `04_marker_dotplot_featureplot.R`, and `cluster9_10_exclusion_suppfig.R` |
| Table EV8 | Top 10 marker genes for neutrophil subclusters 0–10 | `scripts/scRNA/08_neutrophil_marker_dotplot.R` |
| Dataset EV7 | Complete pre-ranked GSEA results for clusters 0–8 | `scripts/gsea/comprehensive/cluster0_GSEA.R`–`cluster8_GSEA.R`; `scripts/gsea/heatmap/02_select_representative_pathways.R` |
