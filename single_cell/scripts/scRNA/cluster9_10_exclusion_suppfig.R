###############################################################################
# Fig. EV4B/C: クラスター9,10 除外根拠の可視化スクリプト
#   Panel A : Cluster 0-10 の nFeature_RNA / nCount_RNA / percent.mt violin plot
#             (中央値を点 + boxplot で重ねる。横3列)  → Cluster 9 除外目的
#   Panel B : Cluster 10 vs Neutrophil clusters(0-8) の marker dot plot
#             (Neutrophil / T / NK / Mast のマーカー)   → Cluster 10 除外目的
#   出力     : PDF(各パネル) + 入力値 CSV
#
# 使い方: 下の【設定】の3項目を自分の環境に合わせて書き換え、全文を実行するだけ。
###############################################################################

## ---- パッケージ -----------------------------------------------------------
# 未インストールなら: install.packages(c("ggplot2","patchwork","dplyr"))
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

## =========================== 【設定】ここだけ変える ========================
# 1) Seurat オブジェクト。RDS から読み込む。
obj <- readRDS(file.path("rds", "07_neutrophil_clustered.rds"))
# 既にメモリ上の obj を使いたい場合は上行をコメントアウトし下行を有効化:
# obj <- get("obj")

# 2) クラスターが入っているメタデータ列名（多くは "seurat_clusters"）
cluster_col <- "seurat_clusters"

# 3) 出力先フォルダ
out_dir <- file.path("Outputs", "FigEV4_cluster9_10_output")
## ==========================================================================

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
DefaultAssay(obj) <- "RNA"

# percent.mt が無ければ計算（マウス: "^mt-"、ヒトなら "^MT-" に変更）
if (!"percent.mt" %in% colnames(obj@meta.data)) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  message("percent.mt をマウスパターン '^mt-' で新規計算しました。")
}

# クラスターを 0..10 の因子として整える
clu <- as.character(obj@meta.data[[cluster_col]])
lv  <- as.character(sort(as.integer(unique(clu))))
obj$.cluster <- factor(clu, levels = lv)

###############################################################################
# Panel A : QC violin (Cluster 0-10)
###############################################################################
qc_metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
stopifnot(all(qc_metrics %in% colnames(obj@meta.data)))

dfA <- data.frame(
  cluster      = obj$.cluster,
  nFeature_RNA = obj@meta.data$nFeature_RNA,
  nCount_RNA   = obj@meta.data$nCount_RNA,
  percent.mt   = obj@meta.data$percent.mt,
  check.names  = FALSE
)

# --- 入力値 CSV(1): 各細胞の QC 値 ---
write.csv(dfA, file.path(out_dir, "FigEV4B_QC_per_cell.csv"), row.names = FALSE)

# --- 入力値 CSV(2): クラスター別サマリー(中央値・平均・四分位) ---
summ <- dfA %>%
  tidyr::pivot_longer(-cluster, names_to = "metric", values_to = "value") %>%
  group_by(metric, cluster) %>%
  summarise(n = n(), mean = mean(value), median = median(value),
            q25 = quantile(value, .25), q75 = quantile(value, .75),
            .groups = "drop") %>%
  arrange(metric, cluster)
write.csv(summ, file.path(out_dir, "FigEV4B_QC_cluster_summary.csv"), row.names = FALSE)

# --- violin 1枚を描く関数(violin + 細いboxplot + 中央値の点) ---
make_vln <- function(metric, ylab) {
  ggplot(dfA, aes(x = cluster, y = .data[[metric]], fill = cluster)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = .2, color = "grey30") +
    geom_boxplot(width = .12, outlier.shape = NA, fill = "white",
                 linewidth = .25, color = "grey20") +
    stat_summary(fun = median, geom = "point", size = 1.6,
                 color = "black", shape = 18) +   # 中央値(◆)
    labs(x = "Cluster", y = ylab, title = metric) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = .5, size = 11))
}

pA1 <- make_vln("nFeature_RNA", "nFeature_RNA")
pA2 <- make_vln("nCount_RNA",   "nCount_RNA")
pA3 <- make_vln("percent.mt",   "percent.mt (%)")

panelA <- pA1 + pA2 + pA3 + plot_layout(ncol = 3)
ggsave(file.path(out_dir, "FigEV4B_QC_violin.pdf"), panelA,
       width = 12, height = 4, useDingbats = FALSE)

###############################################################################
# Panel B : marker dot plot (Cluster 10 vs Neutrophil clusters 0-8)
#   ※ Cluster 9 は別途 Panel A で除外するため B には含めない
###############################################################################

# Marker panel used to assess cluster 10 lineage identity; only genes present in the object are plotted.
marker_list <- list(
  Neutrophil = c("Ly6g","Csf3r","S100a8","S100a9","Retnlg","Ngp","Mpo","Cxcr2","Il1b"),
  `T cell`   = c("Cd3d","Cd3e","Trbc1","Trbc2","Gimap7"),
  NK         = c("Ncr1","Nkg7","Klrd1","Klrb1c","Klri1","Prf1","Gzmb"),
  Mast       = c("Kit","Cpa3","Fcer1a","Cma1","Mcpt4")
)

present <- rownames(obj)
marker_list <- lapply(marker_list, function(g) g[g %in% present])
missing <- setdiff(unlist(list(
  c("Ly6g","Csf3r","S100a8","S100a9","Retnlg","Ngp","Mpo","Cxcr2","Il1b"),
  c("Cd3d","Cd3e","Trbc1","Trbc2","Gimap7"),
  c("Ncr1","Nkg7","Klrd1","Klrb1c","Klri1","Prf1","Gzmb"),
  c("Kit","Cpa3","Fcer1a","Cma1","Mcpt4"))), present)
if (length(missing)) message("オブジェクトに無く除外した遺伝子: ",
                             paste(missing, collapse = ", "))
features <- unlist(marker_list, use.names = FALSE)

# 2群メタデータを作成（Cluster 9 は NA にして除外）
grp <- ifelse(as.character(obj$.cluster) == "10", "Cluster 10",
        ifelse(as.character(obj$.cluster) %in% as.character(0:8),
               "Neutrophil clusters 0-8", NA))
obj$B_group <- factor(grp, levels = c("Neutrophil clusters 0-8", "Cluster 10"))
objB <- subset(obj, cells = colnames(obj)[!is.na(obj$B_group)])
Idents(objB) <- "B_group"

# DotPlot（scale=FALSE: 色を発現量に対応させる。2群だけだとscaled色が±0.707で
#          潰れるのを回避。色は log1p(avg.exp)、点サイズは陽性率 pct.exp）
dp <- DotPlot(objB, features = features, cols = c("lightgrey","firebrick"),
              dot.scale = 6, scale = FALSE) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "firebrick",
                       name = "Average expression\n(log1p)") +
  labs(x = NULL, y = NULL,
       title = "Marker expression: Neutrophil clusters (0-8) vs Cluster 10")

# --- 入力値 CSV(3): dot plot の元数値 + 細胞タイプ注釈 ---
gene2type <- stack(lapply(marker_list, function(g) g))
names(gene2type) <- c("gene", "cell_type")
dpd <- dp$data
dpd$cell_type <- gene2type$cell_type[match(dpd$features.plot, gene2type$gene)]
dpd <- dpd[, c("id","features.plot","cell_type",
               "avg.exp","avg.exp.scaled","pct.exp")]
colnames(dpd)[1:2] <- c("group","gene")
write.csv(dpd, file.path(out_dir, "FigEV4C_marker_dotplot_input.csv"), row.names = FALSE)

ggsave(file.path(out_dir, "FigEV4C_marker_dotplot.pdf"), dp,
       width = max(8, length(features) * .35 + 3), height = 3.5,
       useDingbats = FALSE)

###############################################################################
message("完了。出力先: ", normalizePath(out_dir))
message("  PDF : FigEV4B_QC_violin.pdf / FigEV4C_marker_dotplot.pdf")
message("  CSV : FigEV4B_QC_per_cell.csv / FigEV4B_QC_cluster_summary.csv / FigEV4C_marker_dotplot_input.csv")
###############################################################################
