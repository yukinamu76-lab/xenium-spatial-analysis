## ============================================================
## 03_clustering_umap.R
## Harmony による共有低次元表現を用いたクラスタリング・UMAP
##   → 解析ステップ ① 全細胞の UMAP
## ============================================================
source(file.path("scripts", "scRNA", "00_setup.R"))

clust_dir <- file.path(dir_output, "02_clustering")

merged <- readRDS(file.path(dir_rds, "02_pca.rds"))

## ---- Harmony で sample_id を用いた共有低次元表現を取得 ----
merged <- RunHarmony(merged, group.by.vars = "sample_id", plot_convergence = FALSE)

## ---- ★要設定★ クラスタリングパラメータ ----
n_dims     <- 15    # ★ 02 の Elbow plot を見て決定（例: PC15 あたりから平坦なら 15）★
resolution <- 0.6   # ★ クラスタが粗い/細かすぎる場合は調整（大きいほど細かく分かれる）★

merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:n_dims)
merged <- FindClusters(merged, resolution = resolution)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:n_dims)

## ---- ① UMAP: クラスタ ----
p_umap_cluster <- DimPlot(merged, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("Clusters")
save_plot(p_umap_cluster, "02_umap_clusters.pdf", clust_dir, width = 7, height = 6)

## ---- UMAP: 群別（重ね書き）----
p_umap_group <- DimPlot(merged, reduction = "umap", group.by = "group", cols = group_colors) +
  ggtitle("Group")
save_plot(p_umap_group, "03_umap_group.pdf", clust_dir, width = 7, height = 6)

## ---- UMAP: 群ごとに分割表示 ----
p_umap_split <- DimPlot(merged, reduction = "umap", split.by = "group", label = TRUE, repel = TRUE)
save_plot(p_umap_split, "04_umap_split_by_group.pdf", clust_dir, width = 12, height = 6)

## ---- UMAP: サンプル別 ----
p_umap_sample <- DimPlot(merged, reduction = "umap", group.by = "sample_id")
save_plot(p_umap_sample, "05_umap_sample.pdf", clust_dir, width = 7, height = 6)

saveRDS(merged, file.path(dir_rds, "03_clustered.rds"))
message("03_clustering_umap.R 完了: rds/03_clustered.rds を保存しました")
