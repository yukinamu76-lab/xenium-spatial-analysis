## ============================================================
## 06_neutrophil_subset_pca_elbow.R
## Neutrophil クラスタのみを抽出し、正規化・スケーリング・PCA・Elbow plot をやり直す
## （混入細胞を除いた後は可変遺伝子・主成分の構造が変わるため再計算が必要）
## ============================================================
source(file.path("scripts", "scRNA", "00_setup.R"))

sub_dir <- file.path(dir_output, "05_neutrophil_subset_clustering")

merged <- readRDS(file.path(dir_rds, "03_clustered_annotated.rds"))

## ---- Neutrophil クラスタのみを抽出 ----
neut <- subset(merged, subset = cell_class == "Neutrophil")
message("抽出後の細胞数: ", ncol(neut), " / 全体: ", ncol(merged))

## ---- 正規化・スケーリング・PCA をやり直す ----
neut <- NormalizeData(neut)
neut <- FindVariableFeatures(neut, selection.method = "vst", nfeatures = 2000)
neut <- ScaleData(neut, vars.to.regress = "percent.mt")
neut <- RunPCA(neut, npcs = 30)

## ---- Elbow plot ----
p_elbow <- ElbowPlot(neut, ndims = 30)
save_plot(p_elbow, "01_elbow_plot_neutrophil_only.pdf", sub_dir, width = 6, height = 5)

saveRDS(neut, file.path(dir_rds, "06_neutrophil_pca.rds"))
message("06_neutrophil_subset_pca_elbow.R 完了: rds/06_neutrophil_pca.rds を保存しました")
message(sprintf("%s/05_neutrophil_subset_clustering/01_elbow_plot_neutrophil_only.pdf を確認してください", dir_output))
message("★ 次元数を読み取り、07_neutrophil_clustering_umap.R 内の n_dims を編集してください ★")
