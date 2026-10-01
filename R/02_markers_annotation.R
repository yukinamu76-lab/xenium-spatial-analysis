# =============================================================================
# 02_markers_annotation.R
#
# 目的: クラスターマーカーの同定とセルタイプアノテーション
# 入力: 01_session_save/xen_merged.rds
# 出力: 01_session_save/xen_merged_annotated.rds
#       02_markers_annotation/
#         FindAllMarkers_all.csv, top10_markers.csv, top20_markers.csv
#         perCluster_QC.csv, DotPlot_top3.pdf/png
#         Fig9A_UMAP_annotated.pdf/png, UMAP_clusterID.pdf/png
#         CellType_composition.pdf/png/csv, TableEV10.csv
# 目安: 20〜40分
#
# アノテーションは 00_config.R の annotation_map を書き換えて指定します。
# マーカー検定は sketch 細胞のみのサブセット上で行う (Idents に NA を含めない)。
# 細胞数・割合は全細胞 (cluster_all) ベース。
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})
if (!file.exists("00_config.R")) {
  of <- try(sys.frame(1)$ofile, silent = TRUE)
  if (inherits(of, "try-error") || is.null(of))
    stop("setwd() でスクリプトのあるフォルダへ移動してから実行してください。")
  setwd(dirname(normalizePath(of)))
}
source("00_config.R")
dir.create(DIR_ANNOT, showWarnings = FALSE, recursive = TRUE)

xen_merged <- readRDS(RDS_RAW)
stopifnot(all(c("cluster_all", "sample") %in% colnames(xen_merged@meta.data)),
          SKETCH %in% Assays(xen_merged),
          REDUC_FULL %in% Reductions(xen_merged))
md <- xen_merged@meta.data

# =============================================================================
# 1. クラスター毎の QC 指標 (アノテーション判断の材料)
# =============================================================================
qc_tbl <- md %>%
  group_by(cluster = cluster_all) %>%
  summarise(n = n(),
            n_PBS = sum(sample == "LN_PBS"),
            n_Alb = sum(sample == "LN_Alb"),
            nCount_median   = median(nCount_Xenium),
            nFeature_median = median(nFeature_Xenium),
            pct_nFeature_lt20 = round(mean(nFeature_Xenium < 20) * 100, 2),
            cell_area_median  = if ("cell_area" %in% colnames(md))
              round(median(cell_area, na.rm = TRUE), 2) else NA_real_,
            .groups = "drop") %>%
  arrange(as.numeric(as.character(cluster)))
write.csv(qc_tbl, file.path(DIR_ANNOT, "perCluster_QC.csv"), row.names = FALSE)
cat("\n=== per-cluster QC ===\n"); print(as.data.frame(qc_tbl), row.names = FALSE)
cat(sprintf("overall median: nCount = %.0f, nFeature = %.0f\n",
            median(md$nCount_Xenium), median(md$nFeature_Xenium)))

# =============================================================================
# 2. FindAllMarkers (sketch 細胞のみ)
# =============================================================================
sk_cells <- sketch_cells_of(xen_merged)
stopifnot(!anyNA(md[sk_cells, "seurat_clusters"]))

sk <- subset_cells(xen_merged, sk_cells)
DefaultAssay(sk) <- SKETCH
sk[[SKETCH]] <- JoinLayers(sk[[SKETCH]])
sk$seurat_clusters <- droplevels(factor(as.character(sk$seurat_clusters),
                                        levels = levels(xen_merged$cluster_all)))
Idents(sk) <- "seurat_clusters"
stopifnot(!anyNA(Idents(sk)))
n_clust <- nlevels(sk$seurat_clusters)
message(sprintf("FindAllMarkers: %s cells / %d clusters",
                format(ncol(sk), big.mark = ","), n_clust))

markers <- FindAllMarkers(sk, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, test.use = "wilcox",
                          random.seed = SEED, verbose = FALSE)
write.csv(markers, file.path(DIR_ANNOT, "FindAllMarkers_all.csv"), row.names = FALSE)
markers <- markers %>%
  mutate(cluster = factor(as.character(cluster), levels = levels(xen_merged$cluster_all)))

top_n_tbl <- function(n) {
  markers %>% group_by(cluster) %>%
    slice_max(avg_log2FC, n = n, with_ties = FALSE) %>%
    arrange(cluster, desc(avg_log2FC)) %>% ungroup() %>%
    select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)
}
top10 <- top_n_tbl(10); top20 <- top_n_tbl(20)
write.csv(top10, file.path(DIR_ANNOT, "top10_markers.csv"), row.names = FALSE)
write.csv(top20, file.path(DIR_ANNOT, "top20_markers.csv"), row.names = FALSE)

top3_genes <- top10 %>% group_by(cluster) %>%
  slice_max(avg_log2FC, n = 3, with_ties = FALSE) %>% pull(gene) %>% unique()
save_fig(
  DotPlot(sk, features = top3_genes, scale = FALSE) + RotatedAxis() +
    scale_color_gradient(low = "lightgrey", high = "#D6604D") +
    labs(title = "Top 3 marker genes per cluster", x = NULL, y = "Cluster") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)),
  DIR_ANNOT, "DotPlot_top3",
  width = max(12, length(top3_genes) * 0.35 + 4), height = n_clust * 0.45 + 3, dpi = 300)
rm(sk); invisible(gc(FALSE))

# =============================================================================
# 3. アノテーション付与 (00_config.R の annotation_map)
# =============================================================================
missing_lv <- setdiff(levels(xen_merged$cluster_all), names(annotation_map))
if (length(missing_lv) > 0)
  stop("00_config.R の annotation_map に定義のないクラスター: ",
       paste(missing_lv, collapse = ", "))

xen_merged$annotation <- factor(
  unname(annotation_map[as.character(xen_merged$cluster_all)]), levels = cell_order)
stopifnot(!anyNA(xen_merged$annotation))
xen_merged$is_low_quality <- xen_merged$annotation == LOWQ_LABEL

cat("\n=== cells per annotation ===\n")
print(table(xen_merged$annotation, xen_merged$sample))

saveRDS(xen_merged, RDS_ANNOT)
message("Saved: ", RDS_ANNOT)

# =============================================================================
# 4. Fig.9A — 全細胞UMAP (annotation) / クラスター番号版
# =============================================================================
save_fig(
  DimPlot(xen_merged, reduction = REDUC_FULL, group.by = "annotation",
          cols = annot_colors, label = TRUE, label.size = 2.4, repel = TRUE,
          pt.size = PT_SIZE, raster = TRUE, raster.dpi = RASTER_DPI) +
    scale_color_manual(values = annot_colors, name = "Cell type") +
    labs(title = sprintf("Cell type annotation (n = %s cells)",
                         format(ncol(xen_merged), big.mark = ",")),
         x = "UMAP 1", y = "UMAP 2") +
    theme_pub() + theme(legend.position = "right"),
  DIR_ANNOT, "Fig9A_UMAP_annotated", width = 8, height = 6)

save_fig(
  DimPlot(xen_merged, reduction = REDUC_FULL, group.by = "cluster_all",
          label = TRUE, label.size = 3, repel = TRUE, pt.size = PT_SIZE,
          raster = TRUE, raster.dpi = RASTER_DPI) +
    labs(title = sprintf("Seurat clusters (res = %.1f, dims = %d)", RESOLUTION, N_DIMS),
         x = "UMAP 1", y = "UMAP 2") +
    theme_pub() + theme(legend.position = "right"),
  DIR_ANNOT, "UMAP_clusterID", width = 8, height = 6, dpi = 300)

# =============================================================================
# 5. 細胞組成 (PBS vs Alb)
# =============================================================================
df_bar <- as.data.frame(table(Annotation = xen_merged$annotation,
                              Sample     = xen_merged$sample)) %>%
  group_by(Sample) %>% mutate(pct = Freq / sum(Freq) * 100) %>% ungroup()
write.csv(df_bar %>% pivot_wider(names_from = Sample, values_from = c(Freq, pct)),
          file.path(DIR_ANNOT, "CellType_composition.csv"), row.names = FALSE)

save_fig(
  ggplot(df_bar %>% mutate(Annotation = factor(Annotation, levels = rev(cell_order))),
         aes(Sample, pct, fill = Annotation)) +
    geom_bar(stat = "identity", width = 0.6, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = annot_colors, name = "Cell type",
                      guide = guide_legend(reverse = TRUE)) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.02)),
                       name = "Proportion of cells (%)") +
    scale_x_discrete(name = NULL) + theme_pub() + theme(legend.position = "right"),
  DIR_ANNOT, "CellType_composition", width = 6, height = 6)

# =============================================================================
# 6. Table EV10
# =============================================================================
comp <- md %>% count(cluster = cluster_all, sample, name = "n") %>%
  pivot_wider(names_from = sample, values_from = n, values_fill = 0)
for (s in c("LN_Alb", "LN_PBS")) if (!s %in% colnames(comp)) comp[[s]] <- 0
tot_alb <- sum(comp$LN_Alb); tot_pbs <- sum(comp$LN_PBS); tot_all <- tot_alb + tot_pbs

gene_str <- function(tbl, colname) {
  tbl %>% group_by(cluster) %>%
    summarise(v = paste(gene, collapse = "; "), .groups = "drop") %>%
    mutate(Cluster = as.numeric(as.character(cluster))) %>%
    select(Cluster, !!colname := v)
}

supp18 <- comp %>%
  mutate(Total = LN_Alb + LN_PBS, Cluster = as.numeric(as.character(cluster))) %>%
  arrange(Cluster) %>%
  transmute(
    `Cluster`                    = Cluster,
    `Cell-type annotation`       = unname(annotation_map[as.character(Cluster)]),
    `Alb-treated lymph node, n`  = LN_Alb,
    `PBS-treated lymph node, n`  = LN_PBS,
    `Total cells, n`             = Total,
    `All cells (%)`              = round(Total  / tot_all * 100, 3),
    `Alb-treated lymph node (%)` = round(LN_Alb / tot_alb * 100, 3),
    `PBS-treated lymph node (%)` = round(LN_PBS / tot_pbs * 100, 3)) %>%
  left_join(gene_str(top10, "Top 10 cluster marker genes"), by = "Cluster") %>%
  left_join(gene_str(top20, "Top 20 cluster marker genes"), by = "Cluster")

if (!is.null(annotation_evidence))
  supp18$`Annotation evidence` <- unname(annotation_evidence[as.character(supp18$Cluster)])

write.csv(supp18, file.path(DIR_ANNOT, "TableEV10.csv"),
          row.names = FALSE, na = "")
cat("\n=== Table EV10 ===\n")
print(as.data.frame(supp18 %>% select(1:8)), row.names = FALSE)
cat(sprintf("Total: %s cells (Alb %s / PBS %s)\n",
            format(tot_all, big.mark = ","), format(tot_alb, big.mark = ","),
            format(tot_pbs, big.mark = ",")))

# =============================================================================
# 7. 実行記録
# =============================================================================
write_log(DIR_ANNOT, c(
  paste0("rds in         : ", RDS_RAW),
  paste0("rds out        : ", RDS_ANNOT),
  paste0("cells          : ", ncol(xen_merged), " (PBS ", tot_pbs, " / Alb ", tot_alb, ")"),
  paste0("sketched cells : ", length(sk_cells)),
  paste0("clusters       : ", n_clust),
  paste0("cell types     : ", nlevels(xen_merged$annotation)),
  paste0("marker rows    : ", nrow(markers)),
  paste0(LOWQ_LABEL, " : ", sum(xen_merged$is_low_quality),
         sprintf(" (%.2f%%)", mean(xen_merged$is_low_quality) * 100)),
  "",
  "marker 検定: sketch 細胞のサブセット, wilcox, only.pos, min.pct=0.25, logfc.threshold=0.25",
  "細胞数・割合: 全細胞 (cluster_all) ベース",
  "アノテーション: 00_config.R の annotation_map",
  paste0("低品質細胞の下流での扱い: EXCLUDE_LOW_QUALITY = ", EXCLUDE_LOW_QUALITY)))

message("\n=== 02 done ===")
