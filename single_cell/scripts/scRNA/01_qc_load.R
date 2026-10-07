## ============================================================
## 01_qc_load.R
## Manuscript outputs: sample-level QC values contributing to Table EV6
## サンプル読み込み・QC・フィルタリング
## ============================================================
source(file.path("scripts", "scRNA", "00_setup.R"))

qc_dir <- file.path(dir_output, "01_QC")

## ---- 各サンプルを読み込んで Seurat オブジェクトを作成 ----
seu_list <- list()
for (i in seq_len(nrow(sample_info))) {
  sid  <- sample_info$sample_id[i]
  grp  <- sample_info$group[i]
  path <- sample_info$path[i]

  message("Loading: ", sid, " (", path, ")")

  if (grepl("\\.h5$", path)) {
    counts <- Read10X_h5(path)
  } else {
    counts <- Read10X(data.dir = path)
  }

  so <- CreateSeuratObject(
    counts       = counts,
    project      = sid,
    min.cells    = qc_params$min_cells_per_gene,
    min.features = 0
  )
  so$sample_id <- sid
  so$group     <- grp
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = mito_pattern)

  seu_list[[sid]] <- so
}

## ---- フィルタリング前のバイオリンプロット（サンプルごと）----
merged_raw <- merge(seu_list[[1]], y = seu_list[-1], add.cell.ids = names(seu_list))
merged_raw$sample_id <- factor(merged_raw$sample_id, levels = sample_info$sample_id)

p_qc_before <- VlnPlot(
  merged_raw,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample_id",
  pt.size  = 0,
  ncol     = 3
)
save_plot(p_qc_before, "01_QC_violin_before_filter.pdf", qc_dir, width = 12, height = 5)

p_scatter <- FeatureScatter(merged_raw, feature1 = "nCount_RNA", feature2 = "percent.mt",
                            group.by = "sample_id") +
  FeatureScatter(merged_raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                 group.by = "sample_id")
save_plot(p_scatter, "02_QC_scatter_before_filter.pdf", qc_dir, width = 12, height = 5)

## ---- フィルタリング ----
seu_list_filt <- lapply(seu_list, function(so) {
  subset(
    so,
    subset = nFeature_RNA > qc_params$min_features &
             nFeature_RNA < qc_params$max_features &
             percent.mt   < qc_params$max_percent_mt
  )
})

n_before <- sapply(seu_list, ncol)
n_after  <- sapply(seu_list_filt, ncol)
qc_summary <- tibble(
  sample_id      = names(seu_list),
  group          = as.character(sample_info$group),
  n_cells_before = n_before,
  n_cells_after  = n_after,
  pct_kept       = round(100 * n_after / n_before, 1)
)
write_csv(qc_summary, file.path(qc_dir, "qc_cell_counts.csv"))
print(qc_summary)

## ---- フィルタリング後のバイオリンプロット ----
merged_filt <- merge(seu_list_filt[[1]], y = seu_list_filt[-1], add.cell.ids = names(seu_list_filt))
merged_filt$sample_id <- factor(merged_filt$sample_id, levels = sample_info$sample_id)

p_qc_after <- VlnPlot(
  merged_filt,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample_id",
  pt.size  = 0,
  ncol     = 3
)
save_plot(p_qc_after, "03_QC_violin_after_filter.pdf", qc_dir, width = 12, height = 5)

## ---- 保存（Seurat v5 でレイヤーが分かれている場合は結合しておく）----
if (exists("JoinLayers")) {
  merged_filt <- tryCatch(JoinLayers(merged_filt), error = function(e) merged_filt)
}
saveRDS(merged_filt, file.path(dir_rds, "01_merged_filtered.rds"))
message("01_qc_load.R 完了: rds/01_merged_filtered.rds を保存しました")
