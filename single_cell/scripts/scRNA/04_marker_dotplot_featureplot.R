## ============================================================
## 04_marker_dotplot_featureplot.R
## Manuscript outputs: annotation marker panels contributing to Table EV7 and Fig. EV4A
##   → 解析ステップ ②
##     (a) 免疫細胞タイプ別マーカーの DotPlot（好中球以外の混入クラスタを特定）
##     (b) 好中球マーカーの FeaturePlot（このデータが好中球であることの確認）
##
##   このスクリプトの結果を見て、好中球クラスタ番号を判定し、
##   00_setup.R の neutrophil_clusters に記入してから 05 に進む。
## ============================================================
source(file.path("scripts", "scRNA", "00_setup.R"))

mk_dir <- file.path(dir_output, "03_markers")
merged <- readRDS(file.path(dir_rds, "03_clustered.rds"))

## ---- (a) 免疫細胞タイプ別マーカー（各2〜3遺伝子。必要に応じて変更）----
## T細胞はサブセットまで分けず Tcell 全体のマーカーのみ使用。
immune_cell_markers <- list(
  Neutrophil = c("S100a8", "S100a9", "Ly6g"),
  Monocyte   = c("Ly6c2", "Ccr2", "Csf1r"),
  Macrophage = c("C1qa", "C1qb", "Mrc1"),
  cDC        = c("Clec9a", "Flt3", "Zbtb46"),
  pDC        = c("Siglech", "Bst2", "Tcf4"),
  Basophil   = c("Cpa3", "Gata2", "Prss34"),
  NK         = c("Ncr1", "Nkg7", "Gzmb"),
  Tcell      = c("Cd3d", "Cd3e", "Cd3g")
)

genes_flat <- unlist(immune_cell_markers, use.names = FALSE)
present    <- intersect(genes_flat, rownames(merged))
missing    <- setdiff(genes_flat, rownames(merged))
if (length(missing) > 0) {
  message("データに存在しないため除外された遺伝子: ", paste(missing, collapse = ", "))
}

## 細胞タイプごとにまとまるよう DotPlot の遺伝子順を固定
ordered_genes <- unlist(
  lapply(immune_cell_markers, function(g) intersect(g, present)),
  use.names = FALSE
)

Idents(merged) <- "seurat_clusters"
p_dot <- DotPlot(merged, features = ordered_genes, group.by = "seurat_clusters") +
  RotatedAxis() +
  labs(title = "Immune cell type marker genes across clusters",
       x = NULL, y = "Cluster")
save_plot(p_dot, "01_immune_celltype_dotplot.pdf", mk_dir, width = 12, height = 6)

## ---- (b) 好中球マーカーの FeaturePlot ----
present_neut_markers <- intersect(neutrophil_markers, rownames(merged))
p_neut <- FeaturePlot(merged, features = present_neut_markers, ncol = 3)
save_plot(p_neut, "02_neutrophil_marker_featureplot.pdf", mk_dir, width = 12, height = 8)

message("04_marker_dotplot_featureplot.R 完了")
message(sprintf("%s/03_markers/ の DotPlot・FeaturePlot を確認してください", dir_output))
message("★ 好中球マーカー(S100a8/S100a9/Ly6g/Csf3r/Cxcr2)が高発現のクラスタ番号を判定し、")
message("  00_setup.R の neutrophil_clusters に記入してから 05 に進んでください ★")
