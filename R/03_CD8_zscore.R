# =============================================================================
# 03_CD8_zscore.R
#
# 目的: 全細胞UMAP の出力と、CD8+ T cell の checkpoint / 疲弊マーカー比較
# 入力: 01_session_save/xen_merged_annotated.rds
# 出力: 03_CD8_zscore/
#         LN_UMAP.pdf/png                       (rds からの全細胞UMAP)
#         Fig9B_CD8_Zscore.csv / .pdf           (Fig.9B)
#         Fig9B_CD8_Zscore_barplot.pdf/png
# 目安: 30〜50分
#
# z-score は「全CD8+T細胞 (PBS+Alb 併合) で遺伝子ごとに細胞レベル標準化 →
# 群別に平均」の順で算出する。先に平均発現をとると遺伝子間の発現差が消えるため。
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(patchwork); library(pheatmap); library(grid); library(gridExtra)
})
if (!file.exists("00_config.R")) {
  of <- try(sys.frame(1)$ofile, silent = TRUE)
  if (inherits(of, "try-error") || is.null(of))
    stop("setwd() でスクリプトのあるフォルダへ移動してから実行してください。")
  setwd(dirname(normalizePath(of)))
}
source("00_config.R")
dir.create(DIR_CD8, showWarnings = FALSE, recursive = TRUE)

ln <- readRDS(RDS_ANNOT)
stopifnot(ASSAY %in% Assays(ln), REDUC_FULL %in% Reductions(ln),
          all(c("annotation", "sample") %in% colnames(ln@meta.data)))
if (EXCLUDE_LOW_QUALITY) {
  ln <- subset_cells(ln, colnames(ln)[!ln$is_low_quality])
  message("Low-quality cells を除外: 残り ", ncol(ln), " cells")
}
cat("\n=== annotation ===\n"); print(table(ln$annotation))

# =============================================================================
# 1. 全細胞UMAP
# =============================================================================
cols_use <- annot_colors_for(levels(droplevels(ln$annotation)))

p_anno <- DimPlot(ln, reduction = REDUC_FULL, group.by = "annotation",
                  cols = cols_use, label = TRUE, label.size = 2.4, repel = TRUE,
                  pt.size = PT_SIZE, raster = TRUE, raster.dpi = RASTER_DPI) +
  scale_color_manual(values = cols_use, name = "Cell type") +
  labs(title = sprintf("LN — UMAP (all cells, n = %s)", format(ncol(ln), big.mark = ",")),
       x = "UMAP 1", y = "UMAP 2") +
  theme_pub() + theme(legend.position = "right")

p_smpl <- DimPlot(ln, reduction = REDUC_FULL, group.by = "sample",
                  cols = sample_colors, shuffle = TRUE, pt.size = PT_SIZE,
                  raster = TRUE, raster.dpi = RASTER_DPI) +
  scale_color_manual(values = sample_colors, name = "Sample") +
  labs(title = "Sample", x = "UMAP 1", y = "UMAP 2") + theme_pub()

save_fig(p_anno, DIR_CD8, "LN_UMAP", width = 8, height = 6)
save_fig(p_anno | p_smpl, DIR_CD8, "LN_UMAP_with_sample", width = 15, height = 6)

# =============================================================================
# 2. CD8+ T cell — 細胞レベル z-score → 群平均 (Fig.9B)
# =============================================================================
DefaultAssay(ln) <- ASSAY
ln[[ASSAY]] <- JoinLayers(ln[[ASSAY]])
stopifnot("data" %in% Layers(ln[[ASSAY]]))

meta_l   <- ln@meta.data %>% rownames_to_column("cell_id")
cd8_cells <- meta_l %>% filter(annotation == CD8_LABEL) %>% pull(cell_id)
stopifnot(length(cd8_cells) > 0)
message(sprintf("%s: %s cells", CD8_LABEL, format(length(cd8_cells), big.mark = ",")))

all_genes     <- c(CHECKPOINT_GENES, EXHAUSTION_GENES)
genes_avail   <- intersect(all_genes, rownames(ln[[ASSAY]]))
genes_missing <- setdiff(all_genes, rownames(ln[[ASSAY]]))
if (length(genes_missing) > 0)
  message("Panel に含まれない遺伝子: ", paste(genes_missing, collapse = ", "))
stopifnot(length(genes_avail) > 0)

expr_mat <- as.matrix(
  GetAssayData(ln, assay = ASSAY, layer = "data")[genes_avail, cd8_cells, drop = FALSE])
cell_sample <- setNames(as.character(meta_l$sample), meta_l$cell_id)[colnames(expr_mat)]
stopifnot(!anyNA(cell_sample))

zscore_mat <- t(apply(expr_mat, 1, function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / s
}))
colnames(zscore_mat) <- colnames(expr_mat)

is_pbs <- cell_sample == "LN_PBS"; is_alb <- cell_sample == "LN_Alb"
n_pbs  <- sum(is_pbs); n_alb <- sum(is_alb)
message(sprintf("CD8+ T cells: PBS %s / Alb %s",
                format(n_pbs, big.mark = ","), format(n_alb, big.mark = ",")))

gene_order <- c(intersect(CHECKPOINT_GENES, genes_avail),
                intersect(EXHAUSTION_GENES, genes_avail))
z_pbs <- rowMeans(zscore_mat[gene_order, is_pbs, drop = FALSE], na.rm = TRUE)
z_alb <- rowMeans(zscore_mat[gene_order, is_alb, drop = FALSE], na.rm = TRUE)

p_raw <- sapply(gene_order, function(g)
  suppressWarnings(wilcox.test(zscore_mat[g, is_pbs], zscore_mat[g, is_alb]))$p.value)

result_df <- data.frame(
  Gene              = gene_order,
  Group             = unname(GENE_GROUP[gene_order]),
  PBS_zmean         = round(z_pbs, 4),
  Alb_zmean         = round(z_alb, 4),
  Delta_AlbMinusPBS = round(z_alb - z_pbs, 4),
  PBS_pct_pos       = round(rowMeans(expr_mat[gene_order, is_pbs, drop = FALSE] > 0) * 100, 2),
  Alb_pct_pos       = round(rowMeans(expr_mat[gene_order, is_alb, drop = FALSE] > 0) * 100, 2),
  Wilcoxon_p        = signif(p_raw, 4),
  Wilcoxon_p_BH     = signif(p.adjust(p_raw, method = "BH"), 4),
  n_CD8_PBS         = n_pbs,
  n_CD8_Alb         = n_alb,
  row.names = NULL)

write.csv(result_df, file.path(DIR_CD8, "Fig9B_CD8_Zscore.csv"), row.names = FALSE)
cat("\n=== Fig.9B — CD8+ T cell, mean per-cell z-score ===\n"); print(result_df)

# --- heatmap ------------------------------------------------------------------
z_hm <- as.matrix(result_df[, c("PBS_zmean", "Alb_zmean")])
dimnames(z_hm) <- list(result_df$Gene, c("PBS", "Alb"))
zlim <- max(abs(z_hm), na.rm = TRUE) * 1.1

ph <- pheatmap(
  z_hm,
  color  = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks = seq(-zlim, zlim, length.out = 102),
  cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_row    = data.frame(Group = unname(GENE_GROUP[gene_order]),
                                 row.names = gene_order),
  annotation_colors = list(Group = c(Checkpoint = "#8B1A1A",
                                     `Exhaustion subset` = "#1A3A8B")),
  annotation_names_row = FALSE, border_color = "grey80",
  fontsize = 12, fontsize_row = 11, angle_col = 0,
  display_numbers = TRUE, number_format = "%.3f", fontsize_number = 10,
  silent = TRUE)

pdf(file.path(DIR_CD8, "Fig9B_CD8_Zscore.pdf"), width = 7.5, height = 8)
grid.arrange(
  textGrob("LN CD8+ T: mean per-cell z-score (PBS vs Alb)",
           gp = gpar(fontsize = 14, fontface = "bold")),
  ph$gtable,
  textGrob(sprintf("n = %s CD8+ T cells (PBS %s / Alb %s)",
                   format(n_pbs + n_alb, big.mark = ","),
                   format(n_pbs, big.mark = ","), format(n_alb, big.mark = ",")),
           gp = gpar(fontsize = 9, col = "grey30")),
  ncol = 1, heights = c(0.6, 9.0, 0.4))
dev.off()
message("Saved: Fig9B_CD8_Zscore.pdf")

# --- barplot ------------------------------------------------------------------
save_fig(
  ggplot(result_df %>%
           select(Gene, Group, PBS = PBS_zmean, Alb = Alb_zmean) %>%
           pivot_longer(c(PBS, Alb), names_to = "Sample", values_to = "zmean") %>%
           mutate(Gene = factor(Gene, levels = gene_order),
                  Sample = factor(Sample, levels = c("PBS", "Alb"))),
         aes(Gene, zmean, fill = Sample)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = c(PBS = "#2CA02C", Alb = "#D62728")) +
    facet_grid(~Group, scales = "free_x", space = "free_x") +
    labs(title = "CD8+ T cells — mean per-cell z-score",
         subtitle = sprintf("n = %s (PBS %s / Alb %s)",
                            format(n_pbs + n_alb, big.mark = ","),
                            format(n_pbs, big.mark = ","), format(n_alb, big.mark = ",")),
         x = NULL, y = "mean z-score") +
    theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  DIR_CD8, "Fig9B_CD8_Zscore_barplot", width = 8, height = 5)

# =============================================================================
# 3. 実行記録
# =============================================================================
write_log(DIR_CD8, c(
  paste0("rds            : ", RDS_ANNOT),
  paste0("assay          : ", ASSAY),
  paste0("UMAP reduction : ", REDUC_FULL),
  paste0("cells          : ", ncol(ln)),
  paste0("CD8+ T cells   : ", length(cd8_cells), " (PBS ", n_pbs, " / Alb ", n_alb, ")"),
  paste0("genes used     : ", paste(gene_order, collapse = ", ")),
  paste0("genes missing  : ",
         if (length(genes_missing)) paste(genes_missing, collapse = ", ") else "none"),
  paste0("EXCLUDE_LOW_QUALITY : ", EXCLUDE_LOW_QUALITY),
  "",
  "z-score: 全CD8+T細胞で遺伝子ごとに細胞レベル標準化 → 群別に平均",
  "Wilcoxon 検定: 細胞レベル z-score の PBS vs Alb (BH 補正)"))

message("\n=== 03 done ===")
