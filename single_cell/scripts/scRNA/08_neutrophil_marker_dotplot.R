## ============================================================
## 08_neutrophil_marker_dotplot.R
## Manuscript outputs: top-10 subcluster marker genes contributing to Table EV8
##   → 解析ステップ ④ 好中球マーカーの DotPlot
## Neutrophil クラスタごとの状態マーカー（未成熟/成熟・抑制性・IFN応答・増殖）を可視化
## ============================================================
source(file.path("scripts", "scRNA", "00_setup.R"))

mk_dir <- file.path(dir_output, "06_neutrophil_markers")
neut   <- readRDS(file.path(dir_rds, "07_neutrophil_clustered.rds"))

Idents(neut) <- "seurat_clusters"


## ---- 各クラスタのマーカー遺伝子を網羅的に探索 ----
all_markers <- FindAllMarkers(
  neut,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write_csv(all_markers, file.path(mk_dir, "all_cluster_markers.csv"))

top10 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()

write_csv(top10, file.path(mk_dir, "top10_markers_per_cluster.csv"))

## ---- クラスタ x 好中球サブセット・状態マーカー DotPlot ----
## subset_markers は 00_setup.R で定義（必要に応じてそちらを編集）
present_subset_markers <- intersect(subset_markers, rownames(neut))
missing_subset         <- setdiff(subset_markers, rownames(neut))
if (length(missing_subset) > 0) {
  message("データに存在しないため除外された遺伝子: ", paste(missing_subset, collapse = ", "))
}

p_dot <- DotPlot(neut, features = present_subset_markers, group.by = "seurat_clusters") +
  RotatedAxis() +
  labs(title = "Neutrophil subset/state marker genes across clusters",
       x = NULL, y = "Cluster")
save_plot(p_dot, "01_neutrophil_subset_marker_dotplot.pdf", mk_dir, width = 10, height = 6)

message("08_neutrophil_marker_dotplot.R 完了")
message(sprintf("%s/06_neutrophil_markers/01_neutrophil_subset_marker_dotplot.pdf を確認してください", dir_output))
