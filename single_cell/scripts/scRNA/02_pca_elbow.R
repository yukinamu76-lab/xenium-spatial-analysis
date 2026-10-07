## ============================================================
## 02_pca_elbow.R
## 正規化・スケーリング・PCA・Elbow plot
## （クラスタリングに使う次元数 n_dims を決めるためのチェックポイント）
## ============================================================
source(file.path("scripts", "scRNA", "00_setup.R"))

clust_dir <- file.path(dir_output, "02_clustering")

merged <- readRDS(file.path(dir_rds, "01_merged_filtered.rds"))

## ---- 標準的な前処理 ----
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged, vars.to.regress = "percent.mt")
merged <- RunPCA(merged, npcs = 30)

## ---- Elbow plot ----
p_elbow <- ElbowPlot(merged, ndims = 30)
save_plot(p_elbow, "01_elbow_plot.pdf", clust_dir, width = 6, height = 5)

saveRDS(merged, file.path(dir_rds, "02_pca.rds"))
message("02_pca_elbow.R 完了: rds/02_pca.rds を保存しました")
message(sprintf("%s/02_clustering/01_elbow_plot.pdf を確認してください", dir_output))
message("★ Elbow plot が平坦になる手前の次元数を読み取り、03_clustering_umap.R 内の n_dims を編集してください ★")
