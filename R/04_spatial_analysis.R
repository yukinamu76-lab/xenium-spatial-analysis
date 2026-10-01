# =============================================================================
# 04_spatial_analysis.R
#
# 目的: PMN-MDSC-like cell と exhaustion-associated CD8+ T cell の空間解析
# 入力: 01_session_save/xen_merged_annotated.rds
#       0 rawData/output-*Region_3*/transcripts.parquet (LN_PBS)
#       0 rawData/output-*Region_4*/transcripts.parquet (LN_Alb)
# 出力: 04_spatial/
#         Fig9C_PMNMDSC.csv / _barplot.pdf/png            (Fig.9C)
#         FigEV6_DBSCAN_regions.pdf/png                     (Fig. EV6)
#         Fig9D_NND_Logistic.pdf/png                       (Fig.9D)
#         SpatialMap_overview_*.pdf/png
#         FigEV6_SpatialMap_perRegion.pdf/png           (Fig. EV6)
#         FigEV6_NND_Logistic_perRegion.pdf/png
#         FigEV6_Permutation_test.pdf/png
#         NND_Logistic_summary.csv, Permutation_test_summary.csv
#         NND_CD8_to_PMNMDSC.csv, NND_ViolinInput_PMNtoExhCD8.csv
#         QV_sensitivity.csv, TableEV11.csv
# 目安: 30分〜1時間
#
# 細胞定義 (クラスターに依存しない transcript レベル判定):
#   PMN-MDSC-like : Itgam と Ly6g の transcript がどちらも核内 (全細胞が対象)
#   Exh CD8 T     : Pdcd1 AND (Havcr2 OR Lag3 OR Tigit) が核内 (CD8+ T cells のみ)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(tibble)
  library(arrow); library(dbscan); library(RANN); library(ggplot2); library(patchwork)
})
if (!file.exists("00_config.R")) {
  of <- try(sys.frame(1)$ofile, silent = TRUE)
  if (inherits(of, "try-error") || is.null(of))
    stop("setwd() でスクリプトのあるフォルダへ移動してから実行してください。")
  setwd(dirname(normalizePath(of)))
}
source("00_config.R")
dir.create(DIR_SPATIAL, showWarnings = FALSE, recursive = TRUE)

REG <- region_dirs()
SAMPLES   <- c("LN_PBS", "LN_Alb")
REF_LABEL <- "PMN-MDSC-like"
region_label <- function(lid) paste("Cluster", lid)

# =============================================================================
# 1. 読み込み
# =============================================================================
ln <- readRDS(RDS_ANNOT)
if (EXCLUDE_LOW_QUALITY) {
  ln <- subset_cells(ln, colnames(ln)[!ln$is_low_quality])
  message("Low-quality cells を除外: 残り ", ncol(ln), " cells")
}
meta_l <- ln@meta.data %>% rownames_to_column("cell_id")
stopifnot(all(c("annotation", "sample", "x_centroid", "y_centroid") %in% colnames(meta_l)),
          !anyNA(meta_l$x_centroid))

coords  <- meta_l %>% select(cell_id, x = x_centroid, y = y_centroid, sample, annotation)
cd8_ids <- meta_l %>% filter(annotation == CD8_LABEL) %>% select(cell_id, sample)
all_ids <- meta_l %>% select(cell_id, sample)
n_all   <- sapply(SAMPLES, function(s) sum(meta_l$sample == s))
n_lowq  <- sapply(SAMPLES, function(s)
  sum(meta_l$sample == s & meta_l$annotation == LOWQ_LABEL))

cat("\n=== CD8+ T cells ===\n"); print(table(cd8_ids$sample))

# =============================================================================
# 2. barcode prefix の判定 (cells.parquet の全細胞集合で判定)
# =============================================================================
raw_cell_ids <- function(s) {
  pq  <- file.path(REG[[s]], "cells.parquet")
  csv <- file.path(REG[[s]], "cells.csv.gz")
  if (file.exists(pq)) as.character(arrow::read_parquet(pq, col_select = "cell_id")$cell_id)
  else if (file.exists(csv)) as.character(data.table::fread(csv, select = "cell_id")$cell_id)
  else stop("cells.parquet / cells.csv.gz が見つかりません: ", s)
}

detect_prefix <- function(s) {
  ids  <- all_ids %>% filter(sample == s) %>% pull(cell_id)
  rawi <- raw_cell_ids(s)
  for (p in c(paste0(sub("^LN_", "", s), "_"), paste0(s, "_"), "")) {
    frac <- mean(sub(paste0("^", p), "", ids) %in% rawi)
    message(sprintf("  [%s] prefix '%s' -> %.1f%%", s, p, frac * 100))
    if (frac > 0.9) return(p)
  }
  stop("barcode 一致率が低すぎます: ", s)
}
prefix <- sapply(SAMPLES, detect_prefix, simplify = FALSE)

# =============================================================================
# 3. 核内 transcript の取得 (arrow の遅延評価でフィルタを押し込む)
# =============================================================================
read_nuclear <- function(s) {
  path <- file.path(REG[[s]], "transcripts.parquet")
  message(sprintf("\n[%s] %s", s, path))
  ds  <- arrow::open_dataset(path, format = "parquet")
  nms <- tryCatch(names(ds), error = function(e) ds$schema$names)
  stopifnot(all(c("feature_name", "cell_id", "overlaps_nucleus") %in% nms))
  has_qv <- "qv" %in% nms
  sel <- c("feature_name", "cell_id", "overlaps_nucleus", if (has_qv) "qv")

  pull_rows <- function(do_cast) {
    q <- ds %>% select(all_of(sel))
    if (do_cast) q <- q %>% mutate(feature_name = arrow::cast(feature_name, arrow::string()))
    q %>% filter(overlaps_nucleus == 1, feature_name %in% ALL_NUC_GENES) %>% collect()
  }
  df <- try(pull_rows(FALSE), silent = TRUE)
  if (inherits(df, "try-error")) df <- pull_rows(TRUE)

  df$feature_name <- as.character(df$feature_name)
  df$cell_id      <- as.character(df$cell_id)
  df <- df %>% filter(!is.na(cell_id), cell_id != "UNASSIGNED")
  if (!has_qv) df$qv <- NA_real_
  message(sprintf("  nuclear transcripts: %s rows", format(nrow(df), big.mark = ",")))
  df
}
nuc <- sapply(SAMPLES, read_nuclear, simplify = FALSE)

# 核内カウントの wide 化 → 集団の判定
call_populations <- function(s, qv_min) {
  d <- nuc[[s]]
  if (qv_min > 0 && any(!is.na(d$qv))) d <- d %>% filter(qv >= qv_min)
  wide <- d %>% count(cell_id, feature_name, name = "n") %>%
    pivot_wider(names_from = feature_name, values_from = n, values_fill = 0)
  for (g in ALL_NUC_GENES) if (!g %in% colnames(wide)) wide[[g]] <- 0

  map_ids <- function(raw_pos, seurat_ids)
    seurat_ids[sub(paste0("^", prefix[[s]]), "", seurat_ids) %in% raw_pos]

  list(
    pmn = map_ids(
      wide %>% filter(.data[[PMN_GENES[1]]] >= 1, .data[[PMN_GENES[2]]] >= 1) %>%
        pull(cell_id) %>% as.character(),
      all_ids %>% filter(sample == s) %>% pull(cell_id)),
    exh = map_ids(
      wide %>% filter(.data[[EXH_CORE_GENE]] >= 1,
                      .data[[EXH_OR_GENES[1]]] >= 1 | .data[[EXH_OR_GENES[2]]] >= 1 |
                        .data[[EXH_OR_GENES[3]]] >= 1) %>%
        pull(cell_id) %>% as.character(),
      cd8_ids %>% filter(sample == s) %>% pull(cell_id)))
}

pops    <- sapply(SAMPLES, call_populations, qv_min = QV_MIN, simplify = FALSE)
pmn_ids <- lapply(pops, `[[`, "pmn")
exh_ids <- lapply(pops, `[[`, "exh")

for (s in SAMPLES)
  message(sprintf("%s: %s = %d, Exh CD8 T = %d / %d CD8+T",
                  s, REF_LABEL, length(pmn_ids[[s]]), length(exh_ids[[s]]),
                  sum(cd8_ids$sample == s)))

# QV 閾値の感度解析
sens <- bind_rows(lapply(c(0, 20), function(q) {
  p <- sapply(SAMPLES, call_populations, qv_min = q, simplify = FALSE)
  data.frame(QV_min = q,
             PMN_PBS = length(p$LN_PBS$pmn), PMN_Alb = length(p$LN_Alb$pmn),
             ExhCD8_PBS = length(p$LN_PBS$exh), ExhCD8_Alb = length(p$LN_Alb$exh))
}))
write.csv(sens, file.path(DIR_SPATIAL, "QV_sensitivity.csv"), row.names = FALSE)
cat("\n=== QV sensitivity (本体は QV_MIN =", QV_MIN, ") ===\n")
print(sens, row.names = FALSE)
rm(nuc); invisible(gc(FALSE))

# =============================================================================
# 4. Fig.9C — PMN-MDSC-like の細胞数と割合
# =============================================================================
fig9c <- data.frame(
  Sample      = SAMPLES,
  Total_cells = as.integer(n_all[SAMPLES]),
  PMN_MDSC_n  = sapply(SAMPLES, function(s) length(pmn_ids[[s]])),
  CD8_T_n     = sapply(SAMPLES, function(s) sum(cd8_ids$sample == s)),
  ExhCD8_n    = sapply(SAMPLES, function(s) length(exh_ids[[s]])),
  row.names   = NULL) %>%
  mutate(PMN_MDSC_pct_of_all_cells = round(PMN_MDSC_n / Total_cells * 100, 5),
         PMN_MDSC_per_100k_cells   = round(PMN_MDSC_n / Total_cells * 1e5, 2),
         ExhCD8_pct_of_CD8         = round(ExhCD8_n / CD8_T_n * 100, 3),
         Low_quality_n             = as.integer(n_lowq[Sample]),
         Total_cells_excl_lowQ     = Total_cells - as.integer(n_lowq[Sample]),
         PMN_MDSC_pct_excl_lowQ    = round(PMN_MDSC_n / Total_cells_excl_lowQ * 100, 5))

ft <- fisher.test(matrix(c(fig9c$PMN_MDSC_n, fig9c$Total_cells - fig9c$PMN_MDSC_n), nrow = 2))
fig9c$Fisher_p <- signif(ft$p.value, 4)

write.csv(fig9c, file.path(DIR_SPATIAL, "Fig9C_PMNMDSC.csv"), row.names = FALSE)
cat("\n=== Fig.9C ===\n"); print(as.data.frame(fig9c), row.names = FALSE)
cat(sprintf("Fisher exact test: p = %.4g, OR = %.3f\n", ft$p.value, unname(ft$estimate)))

df9c <- fig9c %>% select(Sample, n = PMN_MDSC_n, pct = PMN_MDSC_pct_of_all_cells) %>%
  mutate(Sample = factor(Sample, levels = SAMPLES))

save_fig(
  (ggplot(df9c, aes(Sample, n, fill = Sample)) +
     geom_col(width = 0.6) + geom_text(aes(label = n), vjust = -0.4, size = 3.5) +
     scale_fill_manual(values = sample_colors) +
     scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
     labs(title = paste(REF_LABEL, "cells"), x = NULL, y = "Number of cells") +
     theme_bw(base_size = 11) + theme(legend.position = "none")) |
    (ggplot(df9c, aes(Sample, pct, fill = Sample)) +
       geom_col(width = 0.6) +
       geom_text(aes(label = sprintf("%.4f%%", pct)), vjust = -0.4, size = 3.2) +
       scale_fill_manual(values = sample_colors) +
       scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
       labs(title = "Proportion of all cells",
            subtitle = sprintf("Fisher exact p = %.3g", ft$p.value),
            x = NULL, y = "% of all cells") +
       theme_bw(base_size = 11) + theme(legend.position = "none")),
  DIR_SPATIAL, "Fig9C_PMNMDSC_barplot", width = 8, height = 4)

# =============================================================================
# 5. DBSCAN による region 分割
# =============================================================================
run_dbscan <- function(s) {
  df <- coords %>% filter(sample == s)
  db <- dbscan::dbscan(as.matrix(df %>% select(x, y)), eps = EPS, minPts = MINPTS)
  df$ln_id <- db$cluster
  cnt   <- table(df$ln_id)
  valid <- setdiff(as.integer(names(cnt[cnt >= MIN_CELLS])), 0L)
  df$ln_id <- ifelse(df$ln_id %in% valid, df$ln_id, 0)
  message(sprintf("[%s] DBSCAN: %d regions (>= %d cells): %s",
                  s, length(valid), MIN_CELLS, paste(valid, collapse = ", ")))
  df
}
coords_by <- sapply(SAMPLES, run_dbscan, simplify = FALSE)

plot_regions <- function(s) {
  ggplot(coords_by[[s]] %>% arrange(ln_id),
         aes(x, y, colour = ifelse(ln_id == 0, "unassigned", region_label(ln_id)))) +
    geom_point(size = 0.15, alpha = 0.5) + coord_equal() +
    labs(title = paste(s, "— DBSCAN regions"), colour = "Region",
         x = "x (um)", y = "y (um)") +
    theme_bw(base_size = 11)
}
save_fig(plot_regions("LN_PBS") | plot_regions("LN_Alb"),
         DIR_SPATIAL, "FigEV6_DBSCAN_regions", width = 14, height = 6)

# =============================================================================
# 6. Spatial map
# =============================================================================
plot_spatial_map <- function(df, s, title) {
  cd8_s <- cd8_ids %>% filter(sample == s) %>% pull(cell_id)
  df <- df %>% mutate(cell_class = case_when(
    cell_id %in% pmn_ids[[s]] ~ REF_LABEL,
    cell_id %in% exh_ids[[s]] ~ "Exh CD8 T",
    cell_id %in% cd8_s        ~ "Other CD8 T",
    TRUE                      ~ "Other"))
  ggplot() +
    geom_point(data = df %>% filter(cell_class == "Other"),
               aes(x, y), colour = "grey85", size = 0.05, alpha = 0.3) +
    geom_point(data = df %>% filter(cell_class == "Other CD8 T"),
               aes(x, y), colour = "#AEC7E8", size = 0.3, alpha = 0.5) +
    geom_point(data = df %>% filter(cell_class == "Exh CD8 T"),
               aes(x, y), colour = "#1F77B4", size = 0.8, alpha = 0.8) +
    geom_point(data = df %>% filter(cell_class == REF_LABEL),
               aes(x, y), colour = "#D62728", size = 2, shape = 17) +
    coord_equal() +
    labs(title = title,
         subtitle = sprintf("%s (red, n=%d) | Exh CD8 T (dark blue, n=%d) | Other CD8 T (light blue)",
                            REF_LABEL, sum(df$cell_class == REF_LABEL),
                            sum(df$cell_class == "Exh CD8 T")),
         x = "x (um)", y = "y (um)") +
    theme_bw(base_size = 11)
}

for (s in SAMPLES)
  save_fig(plot_spatial_map(coords_by[[s]], s, paste(s, "— Spatial map")),
           DIR_SPATIAL, paste0("SpatialMap_overview_", s), width = 9, height = 8)

map_plots <- lapply(
  sort(unique(coords_by[[PRIMARY_SAMPLE]]$ln_id[coords_by[[PRIMARY_SAMPLE]]$ln_id != 0])),
  function(lid) {
    sub <- coords_by[[PRIMARY_SAMPLE]] %>% filter(ln_id == lid)
    plot_spatial_map(sub, PRIMARY_SAMPLE,
                     sprintf("%s — %s (n_%s = %d)", PRIMARY_SAMPLE, region_label(lid),
                             REF_LABEL, sum(sub$cell_id %in% pmn_ids[[PRIMARY_SAMPLE]])))
  })
if (length(map_plots) > 0)
  save_fig(wrap_plots(map_plots, ncol = 2), DIR_SPATIAL, "FigEV6_SpatialMap_perRegion",
           width = 16, height = 7 * ceiling(length(map_plots) / 2), dpi = 300)

# =============================================================================
# 7. NND
#    (A) 各 CD8 T -> 最寄りの PMN-MDSC-like  : logistic regression 用
#    (B) 各 PMN-MDSC-like -> 最寄りの Exh CD8 T : violin plot 用
# =============================================================================
compute_nnd <- function(s) {
  df    <- coords_by[[s]]
  cd8_s <- cd8_ids %>% filter(sample == s) %>% pull(cell_id)
  res_cd8 <- list(); res_ref <- list()
  for (lid in sort(unique(df$ln_id[df$ln_id != 0]))) {
    sub     <- df %>% filter(ln_id == lid)
    pmn_xy  <- sub %>% filter(cell_id %in% pmn_ids[[s]]) %>% select(x, y)
    cd8_sub <- sub %>% filter(cell_id %in% cd8_s)
    exh_xy  <- cd8_sub %>% filter(cell_id %in% exh_ids[[s]]) %>% select(x, y)
    message(sprintf("[%s] %s: %s=%d, CD8 T=%d (Exh=%d)", s, region_label(lid),
                    REF_LABEL, nrow(pmn_xy), nrow(cd8_sub), nrow(exh_xy)))
    if (nrow(pmn_xy) == 0 || nrow(cd8_sub) == 0) next

    nn1 <- RANN::nn2(as.matrix(pmn_xy), as.matrix(cd8_sub %>% select(x, y)), k = 1)
    res_cd8[[region_label(lid)]] <- data.frame(
      ln_id = lid, region = region_label(lid), cell_id = cd8_sub$cell_id,
      dist_um = nn1$nn.dists[, 1],
      cd8_type = ifelse(cd8_sub$cell_id %in% exh_ids[[s]], "Exh CD8 T", "Other CD8 T"))

    if (nrow(exh_xy) > 0) {
      nn2r <- RANN::nn2(as.matrix(exh_xy), as.matrix(pmn_xy), k = 1)
      res_ref[[region_label(lid)]] <- data.frame(
        ln_id = lid, region = region_label(lid), dist_um = nn2r$nn.dists[, 1])
    }
  }
  list(cd8 = if (length(res_cd8)) bind_rows(res_cd8) else
         data.frame(ln_id = integer(0), region = character(0), cell_id = character(0),
                    dist_um = numeric(0), cd8_type = character(0)),
       ref = if (length(res_ref)) bind_rows(res_ref) else
         data.frame(ln_id = integer(0), region = character(0), dist_um = numeric(0)))
}
nnd <- sapply(SAMPLES, compute_nnd, simplify = FALSE)

write.csv(bind_rows(lapply(SAMPLES, function(s)
  if (nrow(nnd[[s]]$cd8)) nnd[[s]]$cd8 %>% mutate(sample = s) else NULL)),
  file.path(DIR_SPATIAL, "NND_CD8_to_PMNMDSC.csv"), row.names = FALSE)

write.csv(bind_rows(lapply(SAMPLES, function(s)
  if (nrow(nnd[[s]]$ref)) nnd[[s]]$ref %>% mutate(sample = s) else NULL)),
  file.path(DIR_SPATIAL, "NND_ViolinInput_PMNtoExhCD8.csv"), row.names = FALSE)

# =============================================================================
# 8. Logistic regression
# =============================================================================
fit_and_plot <- function(df_ln, ref_dist, label) {
  n_exh <- sum(df_ln$cd8_type == "Exh CD8 T")
  if (nrow(df_ln) < 10 || n_exh < 5 || n_exh == nrow(df_ln)) {
    message(sprintf("[%s] skip: n=%d, n_exh=%d", label, nrow(df_ln), n_exh)); return(NULL)
  }
  df_ln <- df_ln %>% mutate(is_exh = as.integer(cd8_type == "Exh CD8 T"),
                            dist_100 = dist_um / 100)
  fit <- glm(is_exh ~ dist_100, data = df_ln, family = binomial)
  b   <- coef(fit)[["dist_100"]]; se <- sqrt(diag(vcov(fit)))[["dist_100"]]
  or  <- exp(b); ci_lo <- exp(b - 1.96 * se); ci_hi <- exp(b + 1.96 * se)
  pv  <- summary(fit)$coefficients["dist_100", "Pr(>|z|)"]

  dr <- data.frame(dist_um = seq(0, quantile(df_ln$dist_um, 0.95), length.out = 300)) %>%
    mutate(dist_100 = dist_um / 100)
  pr <- predict(fit, newdata = dr, type = "response", se.fit = TRUE)
  dr <- dr %>% mutate(prob = pr$fit,
                      lo = pr$fit - 1.96 * pr$se.fit, hi = pr$fit + 1.96 * pr$se.fit)

  obs <- df_ln %>% mutate(bin = cut(dist_um, breaks = 15)) %>% group_by(bin) %>%
    summarise(mid = mean(dist_um), pct = mean(is_exh) * 100, n = n(), .groups = "drop") %>%
    filter(n >= 5)

  p_logit <- ggplot() +
    geom_ribbon(data = dr, aes(dist_um, ymin = lo * 100, ymax = hi * 100),
                fill = "#1F77B4", alpha = 0.2) +
    geom_line(data = dr, aes(dist_um, prob * 100), colour = "#1F77B4", linewidth = 1.1) +
    geom_point(data = obs, aes(mid, pct, size = n), colour = "#1F77B4", alpha = 0.7) +
    labs(title = label,
         subtitle = sprintf("OR/100um=%.3f (95%%CI %.3f-%.3f), p=%.4g | n=%d, n_exh=%d",
                            or, ci_lo, ci_hi, pv, nrow(df_ln), n_exh),
         x = sprintf("Distance to nearest %s (um)", REF_LABEL),
         y = "P(Exh CD8 T) (%)") +
    theme_bw(base_size = 11)

  p_violin <- if (length(ref_dist) >= 1) {
    ggplot(data.frame(dist_um = ref_dist), aes(factor(1), dist_um)) +
      geom_violin(fill = "#2CA02C", alpha = 0.6, scale = "width") +
      geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
      labs(title = paste(label, "— NND:", REF_LABEL, "-> nearest Exh CD8 T"),
           subtitle = sprintf("n = %d (median %.1f um)", length(ref_dist),
                              median(ref_dist)),
           x = NULL, y = "Distance (um)") +
      theme_bw(base_size = 11)
  } else ggplot() + theme_void() + labs(title = paste(label, "— NND"), subtitle = "no data")

  list(plot = p_violin | p_logit, or = or, ci_lo = ci_lo, ci_hi = ci_hi,
       p = pv, n = nrow(df_ln), n_exh = n_exh, label = label)
}

results <- list()
for (s in SAMPLES) {
  d <- nnd[[s]]$cd8
  if (!nrow(d)) next
  for (lid in sort(unique(d$ln_id))) {
    lbl <- sprintf("%s — %s", s, region_label(lid))
    r <- fit_and_plot(d %>% filter(ln_id == lid),
                      nnd[[s]]$ref %>% filter(ln_id == lid) %>% pull(dist_um), lbl)
    if (!is.null(r)) results[[lbl]] <- r
  }
  lbl <- sprintf("%s — Combined", s)
  r <- fit_and_plot(d, nnd[[s]]$ref$dist_um, lbl)
  if (!is.null(r)) results[[lbl]] <- r
}

if (length(results) > 0) {
  save_fig(wrap_plots(lapply(results, `[[`, "plot"), ncol = 1) +
             plot_annotation(title = sprintf("%s -> Exh CD8 T: NND & logistic regression",
                                             REF_LABEL)),
           DIR_SPATIAL, "FigEV6_NND_Logistic_perRegion",
           width = 12, height = 5 * length(results), dpi = 300)

  keys <- grep(paste0("^", PRIMARY_SAMPLE, " — Cluster"), names(results), value = TRUE)
  if (length(keys) > 0)
    save_fig(wrap_plots(lapply(results[keys], `[[`, "plot"), ncol = 1) +
               plot_annotation(title = sprintf("Fig.9D  %s", PRIMARY_SAMPLE)),
             DIR_SPATIAL, "Fig9D_NND_Logistic",
             width = 12, height = 5 * length(keys))

  logit_summary <- bind_rows(lapply(results, function(r)
    data.frame(Label = r$label, OR_per100um = r$or, CI_lo = r$ci_lo, CI_hi = r$ci_hi,
               p_value = r$p, n = r$n, n_exh = r$n_exh)))
  write.csv(logit_summary, file.path(DIR_SPATIAL, "NND_Logistic_summary.csv"),
            row.names = FALSE)
  cat("\n=== logistic regression ===\n")
  print(as.data.frame(logit_summary), row.names = FALSE)
}

# =============================================================================
# 9. Permutation test (Exh ラベルを CD8 T 内でシャッフル)
# =============================================================================
fit_beta <- function(y, x) coef(glm(y ~ x, family = binomial))[["x"]]

run_perm <- function(df, label) {
  df <- df %>% mutate(is_exh = as.integer(cd8_type == "Exh CD8 T"),
                      dist_100 = dist_um / 100)
  n <- nrow(df); n_exh <- sum(df$is_exh)
  if (n_exh < 5 || n_exh == n) {
    message(sprintf("%s: skip permutation (n=%d, n_exh=%d)", label, n, n_exh)); return(NULL)
  }
  obs <- fit_beta(df$is_exh, df$dist_100)
  message(sprintf("%s: n=%d, n_exh=%d, obs_beta=%.4f [%d perms]",
                  label, n, n_exh, obs, N_PERM))
  null_b <- replicate(N_PERM, {
    p <- integer(n); p[sample(n, n_exh)] <- 1L; fit_beta(p, df$dist_100)
  })
  list(label = label, obs = obs, null_b = null_b, n = n, n_exh = n_exh,
       p_two = mean(abs(null_b) >= abs(obs)),
       p_one = if (obs < 0) mean(null_b <= obs) else mean(null_b >= obs))
}

set.seed(SEED)
perms <- list()
for (s in SAMPLES) {
  d <- nnd[[s]]$cd8
  if (!nrow(d)) next
  for (lid in sort(unique(d$ln_id))) {
    lbl <- sprintf("%s — %s", s, region_label(lid))
    r <- run_perm(d %>% filter(ln_id == lid), lbl); if (!is.null(r)) perms[[lbl]] <- r
  }
  lbl <- sprintf("%s — Combined", s)
  r <- run_perm(d, lbl); if (!is.null(r)) perms[[lbl]] <- r
}

if (length(perms) > 0) {
  fmt_p <- function(p) if (p == 0) sprintf("p < %.4f", 1 / N_PERM) else sprintf("p = %.4f", p)
  plots <- lapply(perms, function(r) {
    col <- if (r$p_two < 0.05) "#8B0000" else "grey30"
    ggplot(data.frame(beta = r$null_b), aes(beta)) +
      geom_histogram(bins = 40, fill = "grey70", colour = "white", alpha = 0.9) +
      geom_vline(xintercept = r$obs, colour = col, linewidth = 1.2) +
      annotate("text", x = r$obs, y = Inf, vjust = 1.4, hjust = -0.05,
               label = sprintf("Observed beta=%.4f\nTwo-sided %s\nOne-sided %s",
                               r$obs, fmt_p(r$p_two), fmt_p(r$p_one)),
               size = 3.2, colour = col, fontface = "bold") +
      labs(title = r$label,
           subtitle = sprintf("n=%d, n_Exh=%d | %d permutations", r$n, r$n_exh, N_PERM),
           x = "beta (per 100 um)", y = "Count") +
      theme_bw(base_size = 10)
  })
  save_fig(wrap_plots(plots, ncol = 2) +
             plot_annotation(title = sprintf("Permutation test: %s vs Exh CD8 T", REF_LABEL),
                             subtitle = sprintf("H0: Exh labels random among CD8 T | %d permutations",
                                                N_PERM)),
           DIR_SPATIAL, "FigEV6_Permutation_test",
           width = 12, height = 4 * ceiling(length(plots) / 2), dpi = 300)

  perm_summary <- bind_rows(lapply(perms, function(r)
    data.frame(Label = r$label, obs_beta = r$obs, OR_per100um = exp(r$obs),
               p_two_sided = r$p_two, p_one_sided = r$p_one, n = r$n, n_exh = r$n_exh)))
  write.csv(perm_summary, file.path(DIR_SPATIAL, "Permutation_test_summary.csv"),
            row.names = FALSE)
  cat("\n=== permutation test ===\n"); print(as.data.frame(perm_summary), row.names = FALSE)
}

# =============================================================================
# 10. Table EV11
#     solidity = 50um グリッド占有面積 / 凸包面積 (geometry 判断の目安)
# =============================================================================
hull_area_mm2 <- function(xy) {
  if (nrow(xy) < 3) return(NA_real_)
  h <- chull(xy$x, xy$y); h <- c(h, h[1])
  x <- xy$x[h]; y <- xy$y[h]
  abs(sum(x[-length(x)] * y[-1] - x[-1] * y[-length(y)])) / 2 / 1e6
}

build_region_table <- function(s) {
  df <- coords_by[[s]]
  bind_rows(lapply(sort(unique(df$ln_id[df$ln_id != 0])), function(lid) {
    sub <- df %>% filter(ln_id == lid); xy <- sub %>% select(x, y)
    hull <- hull_area_mm2(xy)
    occ  <- length(unique(paste(floor(xy$x / 50), floor(xy$y / 50)))) * 50^2 / 1e6
    data.frame(
      Sample = s,
      Region = region_label(lid),
      `Area (mm2)`                             = round(hull, 3),
      `Total cells (n)`                        = nrow(sub),
      `PMN-MDSC-like neutrophils (n)`          = sum(sub$cell_id %in% pmn_ids[[s]]),
      `CD8+ T cells (n)`                       = sum(sub$annotation == CD8_LABEL),
      `Exhaustion-associated CD8+ T cells (n)` = sum(sub$cell_id %in% exh_ids[[s]]),
      `Low-quality cells (n)`                  = sum(sub$annotation == LOWQ_LABEL),
      `Occupied area, 50um grid (mm2)`         = round(occ, 3),
      `Solidity (occupied/hull)`               = round(occ / hull, 3),
      `Geometry`                               = NA_character_,
      `Analysis status`                        = NA_character_,
      check.names = FALSE)
  }))
}

supp19_all <- bind_rows(lapply(SAMPLES, build_region_table))
write.csv(supp19_all %>% filter(Sample == PRIMARY_SAMPLE) %>% select(-Sample),
          file.path(DIR_SPATIAL, "TableEV11.csv"), row.names = FALSE, na = "")
write.csv(supp19_all, file.path(DIR_SPATIAL, "TableEV11_allSamples.csv"),
          row.names = FALSE, na = "")
cat("\n=== Table EV11 (", PRIMARY_SAMPLE, ") ===\n", sep = "")
print(as.data.frame(supp19_all %>% filter(Sample == PRIMARY_SAMPLE) %>% select(-Sample)),
      row.names = FALSE)
cat("Geometry / Analysis status は目視判断のため空欄 (Solidity を目安に記入)\n")

# =============================================================================
# 11. 実行記録
# =============================================================================
write_log(DIR_SPATIAL, c(
  paste0("rds             : ", RDS_ANNOT),
  paste0("raw dir         : ", RAW_DIR),
  paste0("primary sample  : ", PRIMARY_SAMPLE),
  paste0("QV_MIN          : ", QV_MIN),
  paste0("DBSCAN          : eps=", EPS, "um, minPts=", MINPTS, ", region >= ", MIN_CELLS, " cells"),
  paste0("permutations    : ", N_PERM),
  paste0("SEED            : ", SEED),
  paste0("prefix          : ", paste(SAMPLES, unlist(prefix), sep = "='", collapse = "' / "), "'"),
  paste0(REF_LABEL, "  : ", paste(sapply(SAMPLES, function(s)
    sprintf("%s %d", s, length(pmn_ids[[s]]))), collapse = " / ")),
  paste0("Exh CD8 T       : ", paste(sapply(SAMPLES, function(s)
    sprintf("%s %d", s, length(exh_ids[[s]]))), collapse = " / ")),
  paste0("EXCLUDE_LOW_QUALITY : ", EXCLUDE_LOW_QUALITY),
  "",
  paste0(REF_LABEL, " : ", paste(PMN_GENES, collapse = " AND "), " の transcript が核内 (全細胞)"),
  paste0("Exh CD8 T      : ", EXH_CORE_GENE, " AND (",
         paste(EXH_OR_GENES, collapse = " OR "), ") が核内 (CD8+ T cells)")))

message("\n=== 04 done. Output: ", DIR_SPATIAL, " ===")
