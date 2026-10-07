## ============================================================
## neutrophil_DEG_full_ranking_exclude9_10.R  （単体実行版）
##   好中球のみ・絞り込みなし DEG（全遺伝子ランキング）の csv 出力
##   【クラスタ 9, 10 を除外して実行する版】
##
##   各クラスタ vs 残りクラスタで FindMarkers を「閾値なし(logfc.threshold=0,
##   min.pct=0)」で実行し、検出された全遺伝子を avg_log2FC 降順でランキングして
##   csv 出力する。この csv を後続の網羅的 GSEA の入力に使う。
##
##   ★9,10 の扱い（下記 2 段階で除外）:
##     (A) ランキングを計算する「対象クラスタ」から 9,10 を外す
##         （＝ cluster9_*.csv / cluster10_*.csv は出力しない）
##     (B) 各クラスタの比較相手である「残り(vs rest)」からも 9,10 を外す
##         （＝ 9,10 の細胞は背景としても使わない）
##     この設定は、既存の neutrophil_cluster_DEGs.R /
##     neutrophil_GSEA_comprehensive_all_clusters.R の除外方針
##     （exclude_cl <- c("9","10")）と同一です。
##
## 【注意】これは有意差で絞り込んだ「DEGリスト」ではなく、閾値なしで検出された
##         全遺伝子のランキングである（p_val_adj 等での絞り込みはしていない）。
##
## 使い方:
##   1. 下の CONFIG の 3 項目（rds_path / out_dir / cluster_col）を環境に合わせて修正
##   2. Rscript neutrophil_DEG_full_ranking_exclude9_10.R  で実行
## ============================================================

library(Seurat)
library(dplyr)
library(readr)

## ---------------- CONFIG（ここだけ環境に合わせて修正）----------------
## 好中球サブセットをクラスタリング済みの Seurat オブジェクト(.rds)
rds_path    <- file.path("rds", "07_neutrophil_clustered.rds")
## 出力先フォルダ（無ければ自動作成）
out_dir     <- file.path("Outputs", "08_neutrophil_DEG")
## クラスタ番号が入っている meta.data の列名
cluster_col <- "seurat_clusters"
## 除外するクラスタ（対象からも「残り(vs rest)」からも外す）
exclude_cl  <- c("9", "10")
## ---------------------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

neut <- readRDS(rds_path)

## Seurat v5 でサンプルごとに layer が分かれている場合は結合しておく
## （これをしないと FindMarkers が layer 単位のままで正しく比較できないことがある）
if ("Assay5" %in% class(neut[[DefaultAssay(neut)]])) {
  neut <- JoinLayers(neut)
}

Idents(neut) <- neut@meta.data[[cluster_col]]

all_clusters <- levels(Idents(neut))
if (is.null(all_clusters)) all_clusters <- sort(unique(as.character(Idents(neut))))

## ランキングを計算する対象クラスタ = 全クラスタから 9,10 を除いたもの
## （このデータでは 0〜8 の 9 クラスタになる）
clusters <- setdiff(all_clusters, exclude_cl)

## ---- 実行前の確認ログ（意図通り 9,10 が除外されているか目視できる）----
cat("検出された全クラスタ  :", paste(all_clusters, collapse = ", "), "\n")
cat("除外クラスタ (9,10)   :", paste(intersect(exclude_cl, all_clusters), collapse = ", "), "\n")
cat("ランキング対象クラスタ:", paste(clusters, collapse = ", "), "\n")

## 念のため：exclude_cl に指定した番号が実際に存在するかを警告
missing_excl <- setdiff(exclude_cl, all_clusters)
if (length(missing_excl) > 0) {
  warning("exclude_cl に指定したが実データに存在しないクラスタ: ",
          paste(missing_excl, collapse = ", "),
          "（クラスタ番号は resolution / n_dims / データで変わります。要確認）")
}

all_rankings <- list()

for (cl in clusters) {
  cat("\n=== Cluster", cl, "(vs rest; 9,10 を除外) ===\n")

  ## 「残り」= 全クラスタから、対象クラスタ自身と exclude_cl(9,10) を除いたもの
  rest_cl <- setdiff(all_clusters, c(cl, exclude_cl))

  markers <- tryCatch(
    FindMarkers(
      neut,
      ident.1         = cl,
      ident.2         = rest_cl,   # ← 9,10 を含まない「残り」
      only.pos        = FALSE,
      min.pct         = 0,         # 閾値なし
      logfc.threshold = 0          # 閾値なし
    ),
    error = function(e) { message("FindMarkers failed for cluster ", cl, ": ", e$message); NULL }
  )

  if (is.null(markers) || nrow(markers) == 0) {
    cat("  結果なし\n"); next
  }

  markers$gene    <- rownames(markers)
  markers$cluster <- cl

  df <- markers %>%
    distinct(gene, .keep_all = TRUE) %>%
    arrange(desc(avg_log2FC)) %>%
    relocate(gene, cluster, .before = everything())

  write_csv(df, file.path(out_dir, paste0("cluster", cl, "_full_gene_ranking.csv")))
  all_rankings[[cl]] <- df
  cat("  ランキング対象遺伝子数:", nrow(df), "\n")
}

## ---- 全クラスタ横断のまとめ（9,10 を含まない）----
combined <- bind_rows(all_rankings)
write_csv(combined, file.path(out_dir, "all_clusters_full_gene_ranking.csv"))

message("\nneutrophil_DEG_full_ranking_exclude9_10.R 完了")
message(sprintf("%s/ に出力しました:", normalizePath(out_dir)))
message("  - cluster*_full_gene_ranking.csv     : クラスタごとの全遺伝子ランキング(閾値なし・9,10除外)")
message("  - all_clusters_full_gene_ranking.csv : 全クラスタ結合版(9,10除外)")
message(sprintf("  ※ 出力対象クラスタ: %s（9,10 は対象・rest の双方から除外済み）",
                paste(clusters, collapse = ", ")))
message("  ※ この csv を後続の網羅的 GSEA の入力に使用する。")
