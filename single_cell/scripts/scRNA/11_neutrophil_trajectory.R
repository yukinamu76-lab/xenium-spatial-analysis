## ============================================================
## 好中球 Trajectory解析 (Slingshot, unsupervised)
##   クラスター9,10は主UMAP塊から孤立しておりtrajectoryに乗らないため、
##   計算前に除外する
## ============================================================

## ---------------- 必要パッケージのインストール ----------------
cran_pkgs <- c("dplyr", "ggplot2", "Seurat", "BiocManager")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

bioc_pkgs <- c("slingshot", "SingleCellExperiment")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)
}

library(Seurat)
library(slingshot)
library(dplyr)
library(ggplot2)

## ---------------- CONFIG ----------------
rds_path         <- file.path("rds", "07_neutrophil_clustered.rds")
cluster_col      <- "seurat_clusters"
group_col        <- "group"
out_dir          <- file.path("Outputs", "09_neutrophil_trajectory")

reduction        <- "umap"
n_dims           <- 2

exclude_clusters <- c("9", "10")   # trajectory計算から除外するクラスター

start_clus       <- NULL           # NULL = 起点を指定せずslingshotに任せる
end_clus         <- NULL           # NULL = 終点も指定しない
## -----------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

so   <- readRDS(rds_path)
meta <- so@meta.data

## ---------------- クラスター除外 ----------------
keep_cells <- rownames(meta)[!as.character(meta[[cluster_col]]) %in% exclude_clusters]
so_sub     <- subset(so, cells = keep_cells)
meta_sub   <- so_sub@meta.data

emb <- Embeddings(so_sub, reduction = reduction)[, 1:n_dims, drop = FALSE]
cl  <- as.character(meta_sub[[cluster_col]])

sds <- slingshot(
  data          = emb,
  clusterLabels = cl,
  start.clus    = start_clus,
  end.clus      = end_clus
)

## ---------------- pseudotime出力 ----------------
pt <- slingPseudotime(sds)
colnames(pt) <- paste0("pseudotime_", colnames(pt))

pt_df <- data.frame(
  cell    = rownames(emb),
  cluster = cl,
  group   = meta_sub[[group_col]],
  pt
)
write.csv(pt_df, file.path(out_dir, "slingshot_pseudotime_per_cell.csv"), row.names = FALSE)

main_pt_col <- grep("^pseudotime_", colnames(pt_df), value = TRUE)[1]
cluster_summary <- pt_df %>%
  group_by(cluster) %>%
  summarise(n = dplyr::n(), mean_pt = mean(.data[[main_pt_col]], na.rm = TRUE),
            median_pt = median(.data[[main_pt_col]], na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_pt)
write.csv(cluster_summary, file.path(out_dir, "slingshot_cluster_mean_pseudotime.csv"), row.names = FALSE)

## ---------------- プロット (全lineage) ----------------
plot_df <- data.frame(x = emb[, 1], y = emb[, 2], cluster = cl, pt = pt_df[[main_pt_col]])

p_cluster <- ggplot(plot_df, aes(x, y, color = cluster)) +
  geom_point(size = 0.6, alpha = 0.7) +
  theme_bw() + labs(x = colnames(emb)[1], y = colnames(emb)[2])

p_pt <- ggplot(plot_df, aes(x, y, color = pt)) +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_color_viridis_c(na.value = "grey80") +
  theme_bw() + labs(x = colnames(emb)[1], y = colnames(emb)[2], color = "pseudotime")

curves <- slingCurves(sds, as.df = TRUE)
p_pt <- p_pt + geom_path(data = curves, aes(x = .data[[names(curves)[1]]], y = .data[[names(curves)[2]]], group = Lineage),
                          color = "black", linewidth = 0.8, inherit.aes = FALSE)

ggsave(file.path(out_dir, "slingshot_cluster_plot.png"), p_cluster, width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "slingshot_pseudotime_plot.png"), p_pt, width = 6.5, height = 5, dpi = 300)

## ---------------- Lineage1のみプロット (PDF, 新規フォルダに出力) ----------------
lineage1_dir <- file.path(out_dir, "lineage1_only")
dir.create(lineage1_dir, showWarnings = FALSE, recursive = TRUE)

# pseudotime_LineageX 列名から Lineage1 に対応する列を特定
pt_lineage1_col <- "pseudotime_Lineage1"
if (!pt_lineage1_col %in% colnames(pt_df)) {
  pt_lineage1_col <- grep("^pseudotime_.*1$", colnames(pt_df), value = TRUE)[1]
}
if (is.na(pt_lineage1_col) || is.null(pt_lineage1_col)) {
  stop("Lineage1に対応するpseudotime列が見つかりません。colnames(pt_df)を確認してください: ",
       paste(colnames(pt_df), collapse = ", "))
}

plot_df_l1 <- data.frame(x = emb[, 1], y = emb[, 2], cluster = cl,
                          pt = pt_df[[pt_lineage1_col]])

# Lineage1のcurveのみ抽出（Lineage列の値は "Lineage1" or 1 などケースに対応）
lineage1_labels <- c("Lineage1", "1", 1)
curves_l1 <- curves[as.character(curves$Lineage) %in% as.character(lineage1_labels), ]

p_pt_l1 <- ggplot(plot_df_l1, aes(x, y, color = pt)) +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_color_viridis_c(na.value = "grey80") +
  geom_path(data = curves_l1,
            aes(x = .data[[names(curves_l1)[1]]], y = .data[[names(curves_l1)[2]]]),
            color = "black", linewidth = 0.8, inherit.aes = FALSE) +
  theme_bw() +
  labs(x = colnames(emb)[1], y = colnames(emb)[2],
       color = "pseudotime (Lineage1)",
       title = "Slingshot Trajectory - Lineage1")

ggsave(file.path(lineage1_dir, "slingshot_lineage1_plot.pdf"), p_pt_l1, width = 6.5, height = 5)

cat("完了:", normalizePath(out_dir), "\n")
cat("Lineage1プロット(PDF)出力:", normalizePath(lineage1_dir), "\n")
