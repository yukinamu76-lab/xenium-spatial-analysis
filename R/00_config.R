# =============================================================================
# Xenium spatial transcriptomics — 4T1 tumor-draining lymph node (PBS vs Alb)
# 00_config.R : 共通設定。01〜04 の各スクリプト冒頭で source() される。
#
# 実行順 (00〜04 を同じフォルダに置いて source() する。
#         source() 経由なら working directory は自動で解決される):
#   source("01_QC_clustering.R")       生データ → QC → クラスタリング → rds
#   source("02_markers_annotation.R")  マーカー同定 → アノテーション → Fig. 9A, Table EV10
#   source("03_CD8_zscore.R")          全細胞UMAP, Fig. 9B
#   source("04_spatial_analysis.R")    Fig. 9C, Fig. 9D, Fig. EV6, Table EV11
#
# 必要パッケージ: Seurat (v5), harmony, arrow, dbscan, RANN, pheatmap,
#                 ggplot2, patchwork, dplyr, tidyr, tibble, gridExtra, data.table
#                 (Seurat の raster 描画では scattermore が利用される場合があります)
#
# 環境に合わせて変更するのは PROJ_ROOT のみ。
# アノテーションは下の annotation_map を書き換えてください。
# =============================================================================

# --- パス --------------------------------------------------------------------
# 公開版: 環境変数 XENIUM_PROJECT_ROOT が設定されていればそれを使用。
# 未設定の場合、R/ ディレクトリから実行していれば1階層上を project root とする。
PROJ_ROOT <- Sys.getenv("XENIUM_PROJECT_ROOT", unset = "")
if (!nzchar(PROJ_ROOT)) {
  wd <- normalizePath(getwd())
  PROJ_ROOT <- if (basename(wd) == "R") dirname(wd) else wd
}

RAW_DIR   <- file.path(PROJ_ROOT, "raw_data")
OUT_ROOT  <- file.path(PROJ_ROOT, "results")

DIR_RDS      <- file.path(OUT_ROOT, "01_session_save")
DIR_QC       <- file.path(OUT_ROOT, "01_QC_clustering")
DIR_ANNOT    <- file.path(OUT_ROOT, "02_markers_annotation")
DIR_CD8      <- file.path(OUT_ROOT, "03_CD8_zscore")
DIR_SPATIAL  <- file.path(OUT_ROOT, "04_spatial")

RDS_RAW   <- file.path(DIR_RDS, "xen_merged.rds")            # 01 の出力
RDS_ANNOT <- file.path(DIR_RDS, "xen_merged_annotated.rds")  # 02 の出力

# --- 解析パラメータ ----------------------------------------------------------
SEED       <- 42

QC_THRESH  <- 3        # nFeature_Xenium > QC_THRESH の細胞を残す
N_SKETCH   <- 50000    # sketch する細胞数 (layer = サンプルごと)
N_DIMS     <- 20       # PCA / Harmony の次元数
RESOLUTION <- 0.5      # FindClusters の resolution

CD8_LABEL  <- "CD8+ T cells"
LOWQ_LABEL <- "Low-quality cells"
EXCLUDE_LOW_QUALITY <- FALSE   # TRUE にすると 03/04 で低品質細胞を除外する

ASSAY        <- "Xenium"
SKETCH       <- "sketch"
REDUC_FULL   <- "umap.full"     # 全細胞UMAP
REDUC_SKETCH <- "umap.sketch"   # sketch 細胞のみのUMAP

# 空間解析 (04)
EPS       <- 50      # DBSCAN eps (um)
MINPTS    <- 100     # DBSCAN minPts
MIN_CELLS <- 5000    # region として採用する最小細胞数
N_PERM    <- 1000    # permutation test の反復数
QV_MIN    <- 0       # transcript の QV 下限 (0 = 無フィルタ)
PRIMARY_SAMPLE <- "LN_PBS"

# CD8 マーカー (03)
CHECKPOINT_GENES <- c("Pdcd1", "Ctla4", "Havcr2", "Lag3", "Tigit")
EXHAUSTION_GENES <- c("Tcf7", "Entpd1", "Cxcl13", "Tox")
GENE_GROUP <- c(
  setNames(rep("Checkpoint",        length(CHECKPOINT_GENES)), CHECKPOINT_GENES),
  setNames(rep("Exhaustion subset", length(EXHAUSTION_GENES)), EXHAUSTION_GENES))

# 空間解析の細胞定義 (04)
PMN_GENES     <- c("Ly6g", "Itgam")                 # 両方が核内 -> PMN-MDSC-like
EXH_CORE_GENE <- "Pdcd1"                            # 核内 AND
EXH_OR_GENES  <- c("Havcr2", "Lag3", "Tigit")       # いずれかが核内 -> Exh CD8 T
ALL_NUC_GENES <- c(PMN_GENES, EXH_CORE_GENE, EXH_OR_GENES)

# --- Region ディレクトリ -----------------------------------------------------
region_dirs <- function() {
  d <- list.dirs(RAW_DIR, recursive = FALSE)
  out <- list(LN_PBS = d[grepl("Region_3", d)][1],
              LN_Alb = d[grepl("Region_4", d)][1])
  stopifnot(!is.na(out$LN_PBS), !is.na(out$LN_Alb))
  out
}

# --- 配色 --------------------------------------------------------------------
sample_colors <- c("LN_PBS" = "#2CA02C", "LN_Alb" = "#D62728")

# 既知のセルタイプ名の固定色。ここに無い名前には ANNOT_PALETTE から自動割り当て。
ANNOT_COLORS_KNOWN <- c(
  "CD8+ T cells"             = "#1F77B4",
  "Tfh cells"                = "#9467BD",
  "Tregs"                    = "#17BECF",
  "NK cells"                 = "#8C564B",
  "Naive/follicular B cells" = "#FF7F0E",
  "GC B cells"               = "#FFBB78",
  "Plasma cells"             = "#BCBD22",
  "Myeloid cells"            = "#D62728",
  "Mature DCs"               = "#98DF8A",
  "cDC1"                     = "#2CA02C",
  "pDCs"                     = "#E7BA52",
  "Stromal cells"            = "#8C8C8C",
  "BECs"                     = "#F7B6D2",
  "LECs"                     = "#E377C2",
  "Pericytes/SMCs"           = "#4D4D4D",
  "Low-quality cells"        = "#DDDDDD")

ANNOT_PALETTE <- c(
  "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1", "#76B7B2",
  "#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC", "#1B9E77", "#D95F02",
  "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")

annot_colors_for <- function(labels) {
  labs <- as.character(unique(labels))
  cols <- setNames(rep(NA_character_, length(labs)), labs)
  hit  <- labs %in% names(ANNOT_COLORS_KNOWN)
  cols[hit] <- unname(ANNOT_COLORS_KNOWN[labs[hit]])
  if (any(is.na(cols))) {
    free <- setdiff(ANNOT_PALETTE, cols[!is.na(cols)])
    if (!length(free)) free <- ANNOT_PALETTE
    cols[is.na(cols)] <- rep(free, length.out = sum(is.na(cols)))
  }
  cols
}

# =============================================================================
# アノテーション  ★ここを書き換えて手動でアノテーションします
#
#   左がクラスター番号、右がセルタイプ名。
#   図の凡例と積み上げ順はこの並び順 (= クラスター番号順) になります。
#   色は上の ANNOT_COLORS_KNOWN にあれば固定色、無ければ自動割り当て。
#
#   判断材料は 02 の出力: top10_markers.csv / top20_markers.csv /
#   DotPlot_top3.pdf / perCluster_QC.csv
#
#   ※ 再クラスタリングでクラスター番号が変わることがあります。
#     マーカーを確認してから書き換えてください。
# =============================================================================
annotation_map <- c(
  "0"  = "CD8+ T cells",
  "1"  = "Naive/follicular B cells",
  "2"  = "Stromal cells",
  "3"  = "Myeloid cells",
  "4"  = "Mature DCs",
  "5"  = "Plasma cells",
  "6"  = "Tregs",
  "7"  = "BECs",
  "8"  = "LECs",
  "9"  = "GC B cells",
  "10" = "Tfh cells",
  "11" = "pDCs",
  "12" = "cDC1",
  "13" = "Low-quality cells",
  "14" = "NK cells",
  "15" = "Pericytes/SMCs"
)

# Table EV10 の 'Annotation evidence' 列。不要なら NULL。
annotation_evidence <- c(
  "0"  = "Cd8a, Cd8b1, Cd3d, Cd3e, Itk, Themis, Lef1, Tcf7, Sell",
  "1"  = "Ighd, Ms4a1, Cd19, Cr2, Cd22, Cd79a, Pax5, Ccr6, Cxcr5",
  "2"  = "Pdgfra, Pdgfrb, Ccl19, Ccl21a, Cxcl13, Col6a3, Col1a2, Vtn",
  "3"  = "Csf1r, Cd68, Mertk, Mafb, Sirpa / S100a8, S100a9, Csf3r, Cxcr2, Cd177",
  "4"  = "Ccr7, Fscn1, Cd80, Flt3, Il15ra, Batf3, Cxcl16, Mmp25",
  "5"  = "Jchain, Mzb1, Sdc1, Prdm1, Igkc, Irf4, Tent5c",
  "6"  = "Foxp3, Il2ra, Ctla4, Ikzf2, Tigit, Tnfrsf4, Izumo1r",
  "7"  = "Pecam1, Cdh5, Sox18, Flt1, Egfl7, Tie1, Plvap, Cd34",
  "8"  = "Prox1, Lyve1, Flt4, Stab2, Tbx1",
  "9"  = "Aicda, S1pr2, Bcl6, Ms4a1, Mki67, Nek2, Pou2af1",
  "10" = "Pdcd1, Il21, Ascl2, Bcl6, Cxcr5, Slamf6, Izumo1r, Cd3d",
  "11" = "Ccr9, Runx2, Irf8, Grm8, Cdh1, Mctp2, Pltp, Clec10a",
  "12" = "Xcr1, Clec9a, Cd207, Tlr11, Itgae, Flt3, Naaa, Id2",
  "13" = "cell-type marker なし。nFeature/nCount/cell_area が全クラスター中最低 (perCluster_QC.csv 参照)",
  "14" = "Ncr1, Gzma, Klrk1, Eomes, Tbx21, Prf1, Il18rap, Il2rb",
  "15" = "Notch3, Pdgfrb, Myh11, Mcam, Des, Itga7, Ednra, Ccn2"
)

cell_order   <- unique(unname(annotation_map))
annot_colors <- annot_colors_for(cell_order)

# --- 作図 --------------------------------------------------------------------
DPI <- 600

# DimPlot の点。raster = TRUE のとき pt.size はラスター画像のピクセル単位で
# 解釈されるため、RASTER_DPI に対して十分大きい値にしないと点が見えなくなる。
RASTER_DPI <- c(2000, 2000)
PT_SIZE    <- 2

# --- 共通関数 ----------------------------------------------------------------
theme_pub <- function(base = 8) {
  ggplot2::theme_classic(base_size = base) +
    ggplot2::theme(
      axis.line     = ggplot2::element_line(linewidth = 0.4),
      axis.ticks    = ggplot2::element_blank(),
      axis.text     = ggplot2::element_blank(),
      legend.key.size = grid::unit(3.2, "mm"),
      legend.text   = ggplot2::element_text(size = base - 0.5),
      legend.title  = ggplot2::element_text(size = base, face = "bold"),
      plot.title    = ggplot2::element_text(size = base + 1, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = base - 0.5, colour = "grey40"))
}

# PDF と PNG を同名で保存
save_fig <- function(plot, dir, name, width, height, dpi = DPI) {
  ggplot2::ggsave(file.path(dir, paste0(name, ".pdf")), plot,
                  width = width, height = height, limitsize = FALSE)
  ggplot2::ggsave(file.path(dir, paste0(name, ".png")), plot,
                  width = width, height = height, dpi = dpi, bg = "white",
                  limitsize = FALSE)
  message("Saved: ", name, ".pdf / .png")
}

# 実行記録
write_log <- function(dir, lines) {
  writeLines(c(paste0("date : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), lines),
             file.path(dir, "run_params.txt"))
  capture.output(sessionInfo(), file = file.path(dir, "sessionInfo.txt"))
}

# sketch assay に含まれる細胞
sketch_cells_of <- function(obj) Cells(obj[[SKETCH]])

# 指定細胞のみのサブセット (FOV が原因で失敗する場合は FOV を外して再試行)
subset_cells <- function(obj, cells) {
  tryCatch(subset(obj, cells = cells), error = function(e) {
    for (im in Images(obj)) obj[[im]] <- NULL
    subset(obj, cells = cells)
  })
}

options(future.globals.maxSize = 8 * 1024^3)
set.seed(SEED)
