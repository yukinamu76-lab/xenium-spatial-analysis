## ============================================================
## 07_neutrophil_clustering_umap.R
## Neutrophil のみのデータで Harmony による共有低次元表現を用いたクラスタリング・UMAP
##   → 解析ステップ ③ 好中球抽出後の UMAP
## （06 の Elbow plot を見てから n_dims を決めて実行）
## ============================================================
source(file.path("scripts", "scRNA", "00_setup.R"))

sub_dir <- file.path(dir_output, "05_neutrophil_subset_clustering")

neut <- readRDS(file.path(dir_rds, "06_neutrophil_pca.rds"))

## ---- Harmony で sample_id を用いた共有低次元表現を取得 ----
neut <- RunHarmony(neut, group.by.vars = "sample_id", plot_convergence = FALSE)

## ---- ★要設定★ クラスタリングパラメータ ----
n_dims     <- 15    # ★ 06 の Elbow plot を見て調整 ★
resolution <- 0.6   # ★ クラスタが粗い/細かすぎる場合は調整 ★

neut <- FindNeighbors(neut, reduction = "harmony", dims = 1:n_dims)
neut <- FindClusters(neut, resolution = resolution)
neut <- RunUMAP(neut, reduction = "harmony", dims = 1:n_dims)

## ---- ③ UMAP: クラスタ ----
p_umap_cluster <- DimPlot(neut, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("Neutrophil clusters")
save_plot(p_umap_cluster, "02_umap_neutrophil_clusters.pdf", sub_dir, width = 7, height = 6)

## ---- UMAP: 群別（重ね書き）----
p_umap_group <- DimPlot(neut, reduction = "umap", group.by = "group", cols = group_colors) +
  ggtitle("Group")
save_plot(p_umap_group, "03_umap_neutrophil_group.pdf", sub_dir, width = 7, height = 6)

## ---- UMAP: 群ごとに分割表示 ----
p_umap_split <- DimPlot(neut, reduction = "umap", split.by = "group", label = TRUE, repel = TRUE)
save_plot(p_umap_split, "04_umap_neutrophil_split_by_group.pdf", sub_dir, width = 12, height = 6)

## ---- UMAP: サンプル別 ----
p_umap_sample <- DimPlot(neut, reduction = "umap", group.by = "sample_id")
save_plot(p_umap_sample, "05_umap_neutrophil_sample.pdf", sub_dir, width = 7, height = 6)

saveRDS(neut, file.path(dir_rds, "07_neutrophil_clustered.rds"))
message("07_neutrophil_clustering_umap.R 完了: rds/07_neutrophil_clustered.rds を保存しました")
