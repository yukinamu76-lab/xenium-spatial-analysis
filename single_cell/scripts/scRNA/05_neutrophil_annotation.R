## ============================================================
## 05_neutrophil_annotation.R
## Manuscript outputs: sample-level neutrophil counts/proportions contributing to Table EV6
## クラスタを Neutrophil / Other に分類し、UMAP で可視化
## （04 の判定に基づく。00_setup.R の neutrophil_clusters を使用）
## ============================================================
source(file.path("scripts", "scRNA", "00_setup.R"))

ann_dir <- file.path(dir_output, "04_neutrophil_annotation")
merged  <- readRDS(file.path(dir_rds, "03_clustered.rds"))

## ---- 好中球クラスタが設定されているか確認 ----
if (is.null(neutrophil_clusters)) {
  stop("neutrophil_clusters が未設定です。04 の結果を見て 00_setup.R に記入してください。")
}

merged$cell_class <- ifelse(
  as.character(merged$seurat_clusters) %in% neutrophil_clusters,
  "Neutrophil", "Other"
)
merged$cell_class <- factor(merged$cell_class, levels = c("Neutrophil", "Other"))

## ---- UMAP: Neutrophil / Other ----
p_class <- DimPlot(merged, reduction = "umap", group.by = "cell_class",
                   cols = c("Neutrophil" = "firebrick", "Other" = "grey70")) +
  ggtitle("Neutrophil vs Other")
save_plot(p_class, "01_umap_neutrophil_vs_other.pdf", ann_dir, width = 7, height = 6)

## ---- サンプルごとの内訳（QC確認用）----
class_summary <- merged@meta.data %>%
  as_tibble() %>%
  count(sample_id, group, cell_class, name = "n_cells") %>%
  group_by(sample_id) %>%
  mutate(total_cells = sum(n_cells), pct = round(100 * n_cells / total_cells, 1)) %>%
  ungroup()
write_csv(class_summary, file.path(ann_dir, "cell_class_summary.csv"))
print(class_summary)

saveRDS(merged, file.path(dir_rds, "03_clustered_annotated.rds"))
message("05_neutrophil_annotation.R 完了: rds/03_clustered_annotated.rds を保存しました")
