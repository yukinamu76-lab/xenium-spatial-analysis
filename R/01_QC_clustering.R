# =============================================================================
# 01_QC_clustering.R
#
# 目的: Xenium 生データ → QC → 統合 → クラスタリング → 全細胞UMAP → rds 保存
# 入力: raw_data/<directory containing Region_3> (LN_PBS) / <directory containing Region_4> (LN_Alb)
# 出力: results/01_session_save/xen_merged.rds
#       results/01_QC_clustering/ (QC図, Elbow, UMAP, run_params.txt, sessionInfo.txt)
# 目安: 1〜2時間
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(arrow)
  library(ggplot2); library(patchwork); library(dplyr)
})
if (!file.exists("00_config.R")) {
  of <- try(sys.frame(1)$ofile, silent = TRUE)
  if (inherits(of, "try-error") || is.null(of))
    stop("setwd() でスクリプトのあるフォルダへ移動してから実行してください。")
  setwd(dirname(normalizePath(of)))
}
source("00_config.R")
dir.create(DIR_RDS, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_QC,  showWarnings = FALSE, recursive = TRUE)

REG <- region_dirs()
message("Region_3 (LN_PBS): ", REG$LN_PBS)
message("Region_4 (LN_Alb): ", REG$LN_Alb)

# =============================================================================
# 1. 読み込み
#    LoadXenium は centroid を meta.data に入れないため cells.parquet から付与する
#    (下流の空間解析が meta.data$x_centroid / y_centroid を使用)
# =============================================================================
read_cells_table <- function(region_dir) {
  pq  <- file.path(region_dir, "cells.parquet")
  csv <- file.path(region_dir, "cells.csv.gz")
  if (file.exists(pq)) {
    as.data.frame(arrow::read_parquet(
      pq, col_select = c("cell_id", "x_centroid", "y_centroid",
                         "cell_area", "nucleus_area")))
  } else if (file.exists(csv)) {
    as.data.frame(data.table::fread(csv))
  } else stop("cells.parquet / cells.csv.gz が見つかりません: ", region_dir)
}

load_region <- function(region_dir, fov_name, ident) {
  message("\n=== Loading ", ident, " ===")
  obj <- LoadXenium(region_dir, fov = fov_name)
  obj$orig.ident <- ident

  tbl <- read_cells_table(region_dir)
  tbl$cell_id <- as.character(tbl$cell_id)
  idx <- match(colnames(obj), tbl$cell_id)
  if (anyNA(idx)) stop("cells テーブルと barcode が一致しません: ", ident)

  obj$x_centroid <- tbl$x_centroid[idx]
  obj$y_centroid <- tbl$y_centroid[idx]
  if ("cell_area"    %in% colnames(tbl)) obj$cell_area    <- tbl$cell_area[idx]
  if ("nucleus_area" %in% colnames(tbl)) obj$nucleus_area <- tbl$nucleus_area[idx]

  message(sprintf("  cells = %d, genes = %d", ncol(obj), nrow(obj)))
  obj
}

xen_pbs <- load_region(REG$LN_PBS, "PBS", "PBS")
xen_alb <- load_region(REG$LN_Alb, "Alb", "Alb")

# =============================================================================
# 2. QC (フィルタ前)
# =============================================================================
df_qc <- bind_rows(
  data.frame(sample = "LN_PBS", nCount = xen_pbs$nCount_Xenium,
             nFeature = xen_pbs$nFeature_Xenium),
  data.frame(sample = "LN_Alb", nCount = xen_alb$nCount_Xenium,
             nFeature = xen_alb$nFeature_Xenium)) %>%
  mutate(sample = factor(sample, levels = names(sample_colors)))

qc_sum <- df_qc %>% group_by(sample) %>%
  summarise(n_total = n(),
            n_pass       = sum(nFeature > QC_THRESH),
            n_removed    = sum(nFeature <= QC_THRESH),
            pct_removed  = round(mean(nFeature <= QC_THRESH) * 100, 3),
            median_nFeature = median(nFeature),
            median_nCount   = median(nCount),
            .groups = "drop")
write.csv(qc_sum, file.path(DIR_QC, "QC_summary.csv"), row.names = FALSE)
cat("\n=== QC summary ===\n"); print(as.data.frame(qc_sum), row.names = FALSE)

vln <- function(y, ttl, ylab) {
  ggplot(df_qc, aes(sample, .data[[y]], fill = sample)) +
    geom_violin(scale = "width", alpha = 0.8) +
    geom_boxplot(width = 0.06, outlier.shape = NA, fill = "white") +
    scale_fill_manual(values = sample_colors) +
    scale_y_continuous(limits = c(0, quantile(df_qc[[y]], 0.99))) +
    labs(title = ttl, x = NULL, y = ylab) +
    theme_bw(base_size = 12) + theme(legend.position = "none")
}

p_hist <- ggplot(df_qc %>% filter(nFeature <= 50), aes(nFeature, fill = sample)) +
  geom_histogram(binwidth = 1, alpha = 0.85, colour = "white", position = "identity") +
  geom_vline(xintercept = QC_THRESH + 0.5, colour = "red", linetype = "dashed") +
  annotate("text", x = QC_THRESH + 2, y = Inf, label = sprintf("cutoff > %d", QC_THRESH),
           colour = "red", hjust = 0, vjust = 1.5, size = 3.5) +
  scale_fill_manual(values = sample_colors) +
  facet_wrap(~sample, ncol = 1, scales = "free_y") +
  labs(title = "nFeature_Xenium (<= 50)", x = "Number of genes", y = "Cell count") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

save_fig(
  (vln("nFeature", "nFeature_Xenium", "Number of genes") +
     geom_hline(yintercept = QC_THRESH, linetype = "dashed", colour = "red") |
     vln("nCount", "nCount_Xenium", "Number of transcripts") | p_hist) +
    plot_annotation(title = "Xenium QC — 4T1 tumor-draining lymph node",
                    subtitle = sprintf("QC: nFeature_Xenium > %d", QC_THRESH)),
  DIR_QC, "QC", width = 14, height = 6, dpi = 300)

# =============================================================================
# 3. フィルタ + 統合
# =============================================================================
xen_pbs <- subset(xen_pbs, subset = nFeature_Xenium > QC_THRESH)
xen_alb <- subset(xen_alb, subset = nFeature_Xenium > QC_THRESH)

xen_merged <- merge(xen_pbs, y = xen_alb,
                    add.cell.ids = c("PBS", "Alb"), merge.data = FALSE)
xen_merged$sample <- factor(
  dplyr::recode(as.character(xen_merged$orig.ident), "PBS" = "LN_PBS", "Alb" = "LN_Alb"),
  levels = names(sample_colors))

message(sprintf("\nMerged: %d cells (PBS %d / Alb %d)",
                ncol(xen_merged), ncol(xen_pbs), ncol(xen_alb)))
rm(xen_pbs, xen_alb); invisible(gc(FALSE))

stopifnot(all(c("x_centroid", "y_centroid") %in% colnames(xen_merged@meta.data)),
          !anyNA(xen_merged$x_centroid))

# =============================================================================
# 4. 正規化 → sketch → PCA
#    JoinLayers はしない (ProjectIntegration が layer 分割を必要とする)
# =============================================================================
DefaultAssay(xen_merged) <- ASSAY
xen_merged <- NormalizeData(xen_merged, normalization.method = "LogNormalize", verbose = FALSE)
xen_merged <- FindVariableFeatures(xen_merged, verbose = FALSE)

set.seed(SEED)
xen_merged <- SketchData(xen_merged, ncells = N_SKETCH, method = "LeverageScore",
                        sketched.assay = SKETCH, seed = SEED)

DefaultAssay(xen_merged) <- SKETCH
xen_merged <- FindVariableFeatures(xen_merged, verbose = FALSE)
xen_merged <- ScaleData(xen_merged, verbose = FALSE)
xen_merged <- RunPCA(xen_merged, npcs = 50, seed.use = SEED, verbose = FALSE)

save_fig(
  ElbowPlot(xen_merged, ndims = 50) +
    geom_vline(xintercept = N_DIMS, colour = "red", linetype = "dashed") +
    annotate("text", x = N_DIMS + 1, y = Inf, label = sprintf("dims = %d", N_DIMS),
             colour = "red", hjust = 0, vjust = 1.5, size = 3.5) +
    labs(title = "Elbow plot (sketch assay)") + theme_bw(base_size = 12),
  DIR_QC, "ElbowPlot", width = 8, height = 5, dpi = 300)

# =============================================================================
# 5. Harmony 統合 → UMAP → クラスタリング (sketch)
# =============================================================================
xen_merged <- IntegrateLayers(xen_merged, method = HarmonyIntegration,
                              orig.reduction = "pca", new.reduction = "harmony",
                              group.by = "orig.ident", verbose = FALSE)

xen_merged <- RunUMAP(xen_merged, reduction = "harmony", dims = 1:N_DIMS,
                      reduction.name = REDUC_SKETCH, reduction.key = "UMAPsketch_",
                      return.model = TRUE, seed.use = SEED, verbose = FALSE)
xen_merged <- FindNeighbors(xen_merged, reduction = "harmony", dims = 1:N_DIMS, verbose = FALSE)
xen_merged <- FindClusters(xen_merged, resolution = RESOLUTION,
                           random.seed = SEED, verbose = FALSE)
n_clust_sketch <- nlevels(xen_merged$seurat_clusters)
message(sprintf("Clusters on sketch (res = %.1f, dims = %d): %d",
                RESOLUTION, N_DIMS, n_clust_sketch))

# =============================================================================
# 6. 全細胞へ投影 → 全細胞UMAP
# =============================================================================
xen_merged <- ProjectIntegration(xen_merged, sketched.assay = SKETCH, assay = ASSAY,
                                 reduction = "harmony", seed = SEED)
xen_merged <- ProjectData(xen_merged, sketched.assay = SKETCH, assay = ASSAY,
                          sketched.reduction = "harmony.full",
                          full.reduction     = "harmony.full",
                          dims = 1:N_DIMS,
                          refdata = list(cluster_full = "seurat_clusters"))
DefaultAssay(xen_merged) <- ASSAY

xen_merged$cluster_all <- ifelse(is.na(as.character(xen_merged$seurat_clusters)),
                                 as.character(xen_merged$cluster_full),
                                 as.character(xen_merged$seurat_clusters))
stopifnot(!anyNA(xen_merged$cluster_all))
xen_merged$cluster_all <- factor(
  xen_merged$cluster_all,
  levels = as.character(sort(as.numeric(unique(xen_merged$cluster_all)))))
message(sprintf("Clusters (all cells): %d", nlevels(xen_merged$cluster_all)))
print(table(xen_merged$cluster_all, xen_merged$sample))

xen_merged <- RunUMAP(xen_merged, reduction = "harmony.full", dims = 1:N_DIMS,
                      reduction.name = REDUC_FULL, reduction.key = "UMAPfull_",
                      seed.use = SEED, verbose = FALSE)

saveRDS(xen_merged, RDS_RAW)
message("Saved: ", RDS_RAW)

# =============================================================================
# 7. UMAP 図
# =============================================================================
save_fig(
  DimPlot(xen_merged, reduction = REDUC_FULL, group.by = "cluster_all",
          label = TRUE, label.size = 3, repel = TRUE, pt.size = PT_SIZE,
          raster = TRUE, raster.dpi = RASTER_DPI) +
    labs(title = sprintf("Seurat clusters (res = %.1f, dims = %d, n = %s cells)",
                         RESOLUTION, N_DIMS, format(ncol(xen_merged), big.mark = ",")),
         x = "UMAP 1", y = "UMAP 2") +
    theme_pub() + theme(legend.position = "right"),
  DIR_QC, "UMAP_cluster", width = 8, height = 6)

save_fig(
  DimPlot(xen_merged, reduction = REDUC_FULL, group.by = "sample",
          cols = sample_colors, shuffle = TRUE, pt.size = PT_SIZE,
          raster = TRUE, raster.dpi = RASTER_DPI) +
    scale_color_manual(values = sample_colors, name = "Sample") +
    labs(title = "Sample", x = "UMAP 1", y = "UMAP 2") + theme_pub(),
  DIR_QC, "UMAP_sample", width = 7, height = 6)

# =============================================================================
# 8. 実行記録
# =============================================================================
write_log(DIR_QC, c(
  paste0("region PBS/Alb   : ", basename(REG$LN_PBS), " / ", basename(REG$LN_Alb)),
  paste0("panel genes      : ", nrow(xen_merged[[ASSAY]])),
  paste0("QC               : nFeature_Xenium > ", QC_THRESH),
  paste0("N_SKETCH         : ", N_SKETCH, " (layer = sample ごと)"),
  paste0("sketched cells   : ", length(sketch_cells_of(xen_merged))),
  paste0("N_DIMS           : ", N_DIMS),
  paste0("RESOLUTION       : ", RESOLUTION),
  paste0("SEED             : ", SEED),
  paste0("cells after QC   : ", ncol(xen_merged),
         " (PBS ", sum(xen_merged$sample == "LN_PBS"),
         " / Alb ", sum(xen_merged$sample == "LN_Alb"), ")"),
  paste0("clusters sketch  : ", n_clust_sketch),
  paste0("clusters all     : ", nlevels(xen_merged$cluster_all)),
  paste0("reductions       : ", paste(Reductions(xen_merged), collapse = ", ")),
  paste0("rds              : ", RDS_RAW)))

message("\n=== 01 done ===")
