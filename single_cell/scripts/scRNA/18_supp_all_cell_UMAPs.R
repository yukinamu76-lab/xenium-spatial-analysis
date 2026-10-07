## ============================================================
## 18_supp_all_cell_UMAPs.R
## Fig. EV4A: 全細胞UMAP（アノテーション前 / アノテーション後）
## 使用データ: 03_clustering_umap.R の全細胞クラスタリング結果 (rds/03_clustered.rds)
##   (好中球サブセット後の07番オブジェクトではなく、全細胞オブジェクトを使用)
## ============================================================
source(file.path("scripts", "scRNA", "00_setup.R"))

out_dir <- file.path(dir_output, "03_all_UMAP")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

merged <- readRDS(file.path(dir_rds, "03_clustered.rds"))

## ---- (1) アノテーション前: クラスタ番号のみのUMAP ----
p_preanno <- DimPlot(merged, reduction = "umap", group.by = "seurat_clusters",
                      label = TRUE, repel = TRUE) +
  labs(title = NULL)
save_plot(p_preanno, "FigEV4A_all_cells_UMAP_preannotation.pdf", out_dir, width = 7, height = 6)

## ---- (2) アノテーション後: 好中球以外は「Other」ではなく各cell typeで表示 ----
## クラスタ→細胞タイプの対応 (04_immune_marker_dotplot.R のマーカーDotPlotに基づく判定。
## 00_setup.R の neutrophil_clusters / コメントと同一の対応関係)
cluster_celltype_map <- setNames(rep("Neutrophil", length(neutrophil_clusters)), neutrophil_clusters)
cluster_celltype_map["10"] <- "Monocyte"
cluster_celltype_map["11"] <- "NK"
cluster_celltype_map["13"] <- "pDC"
cluster_celltype_map["15"] <- "Basophil"

## unname()で外す: 名前付きベクトルのままだとSeuratが「クラスタ番号」を
## セルバーコードとして解釈しようとしエラーになるため、位置対応の代入にする
merged$cell_type <- unname(cluster_celltype_map[as.character(merged$seurat_clusters)])
merged$cell_type <- factor(merged$cell_type,
                            levels = c("Neutrophil", "Monocyte", "NK", "pDC", "Basophil"))

celltype_cols <- c(Neutrophil = "#F0C808", Monocyte = "#1B9E77",
                    NK = "#D95F02", pDC = "#7570B3", Basophil = "#E7298A")

p_annotated <- DimPlot(merged, reduction = "umap", group.by = "cell_type",
                        label = TRUE, repel = TRUE, cols = celltype_cols) +
  labs(title = NULL)
save_plot(p_annotated, "FigEV4A_all_cells_UMAP_annotated_celltype.pdf", out_dir, width = 7.5, height = 6)

## ---- クラスタ→細胞タイプ対応表（annotation reference CSV）----
cluster_table <- tibble(
  cluster = names(cluster_celltype_map),
  cell_type = unname(cluster_celltype_map)
) %>% arrange(as.integer(cluster))
write_csv(cluster_table, file.path(out_dir, "cluster_celltype_table.csv"))

message("18_supp_all_cell_UMAPs.R 完了")
message(sprintf("%s を確認してください", out_dir))
