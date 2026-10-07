## ============================================================
## 00_setup.R
## 好中球 scRNA-seq 解析：共通設定
## （全スクリプトの冒頭で source(file.path("scripts", "scRNA", "00_setup.R")) される）
## ============================================================
##
## 使い方:
##   1. RStudio で本フォルダ（scripts/ の親）をプロジェクトルートとして開く。
##   2. 下記 CONFIG の ★要編集★ 箇所（プロジェクトのパス・サンプル定義・群設定）を
##      自分のデータに合わせて書き換える。
##   3. scripts/ を 01 → 02 →(Elbow確認)→ 03 → 04 → 05 → 06 →(Elbow確認)→
##      07 → 08 → 09 → 10 → 11 の順に実行する。
##   4. 各スクリプトは <project_root>/Outputs 以下に番号付きフォルダを作り、
##      図(pdf)・表(csv)を出力する。中間オブジェクトは <project_root>/rds に保存される。
##
## 出力形式: 図はすべて PDF、表はすべて CSV に統一している。
##
## 解析の全体像（①〜⑥は解析ステップ、括弧内は担当スクリプト）:
##   ① 全細胞の UMAP                          → 01, 02, 03
##   ② マーカー DotPlot + 好中球 FeaturePlot    → 04
##   ③ 好中球抽出後の UMAP                     → 05, 06, 07
##   ④ 好中球マーカー DotPlot                  → 08
##   （+）好中球クラスタの群間 組成比較 csv       → 09（③の後に実行）
##   ⑤ 好中球のみ・絞り込みなし DEG の csv 出力  → 10
##   ⑥ 好中球のみ・起点/終点指定なし Trajectory  → 11
## ------------------------------------------------------------

## ---- パッケージ ----
required_cran <- c("Seurat", "dplyr", "tidyr", "ggplot2", "patchwork",
                   "readr", "tibble", "RColorBrewer", "harmony", "scales")
missing_cran <- setdiff(required_cran, rownames(installed.packages()))
if (length(missing_cran) > 0) install.packages(missing_cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

## 09（組成比較）で使用。レプリケートを考慮した検定 propeller に必要。
## 入らない場合は 09 が Wilcoxon/t 検定のみにフォールバックする。
if (!requireNamespace("speckle", quietly = TRUE)) {
  tryCatch(BiocManager::install("speckle", update = FALSE, ask = FALSE),
           error = function(e) message("speckle のインストールに失敗（09 は Wilcoxon/t 検定のみ使用）"))
}

## 11（Trajectory）で使用。
bioc_pkgs <- c("slingshot", "SingleCellExperiment")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    tryCatch(BiocManager::install(pkg, update = FALSE, ask = FALSE),
             error = function(e) message(pkg, " のインストールに失敗しました（11実行時に再確認）"))
  }
}

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)
library(tibble)
library(RColorBrewer)
library(harmony)

set.seed(1234)   # 乱数固定（再現性のため変更しないことを推奨）

## ============================================================
## ★★★ CONFIG（ここを自分のデータに合わせて編集）★★★
## ============================================================

## ---- ★要編集★ プロジェクトのルート ----
## RStudio Project として開いていれば getwd() がプロジェクトルートになる。
## 直接パスを指定する場合は下行を書き換える。
proj_dir  <- getwd()
## proj_dir <- "/path/to/your/project"   # ← 絶対パス指定したい場合はこちら

## ---- ★要編集★ 入力データの置き場所 ----
## 各サンプルの 10x フィルタ済みマトリクス（barcodes/features/matrix の3ファイル形式、
## または .h5）を置いたフォルダの親ディレクトリ。
data_root <- file.path(proj_dir, "Rawdata")

## ---- ★要編集★ サンプル定義 ----
## sample_id : 任意の一意なサンプル名
## group     : 比較したい群ラベル（2群を想定。名前は任意）
## path      : そのサンプルの10xマトリクスフォルダ（または.h5ファイル）へのパス
sample_info <- tribble(
  ~sample_id, ~group, ~path,
  "PBS-1",   "PBS", file.path(data_root, "PBS_1_sample_filtered_feature_bc_matrix.h5"),
  "PBS-2",   "PBS", file.path(data_root, "PBS_2_sample_filtered_feature_bc_matrix.h5"),
  "Alb-1",   "Alb", file.path(data_root, "Alb_1_sample_filtered_feature_bc_matrix.h5"),
  "Alb-2",   "Alb", file.path(data_root, "Alb_2_sample_filtered_feature_bc_matrix.h5")
)

## ---- ★要編集★ 群の表示順と色 ----
## group_levels の先頭を基準群（対照群）にする。sample_info$group と一致させること。
group_levels <- c("PBS", "Alb")
group_colors <- c("PBS" = "grey60", "Alb" = "firebrick")

sample_info$group <- factor(sample_info$group, levels = group_levels)

## ---- ★要編集★ 生物種（ミトコンドリア遺伝子のプレフィックス）----
## マウス: "^mt-"  /  ヒト: "^MT-"
mito_pattern <- "^mt-"

## ---- ★要編集の可能性あり★ QC閾値 ----
qc_params <- list(
  min_features       = 200,    # この数未満の遺伝子数の細胞を除外
  max_features       = 6000,   # この数を超える遺伝子数の細胞を除外（ダブレット対策）
  max_percent_mt     = 15,     # ミトコンドリア遺伝子割合(%)の上限
  min_cells_per_gene = 3       # この数未満の細胞でしか検出されない遺伝子を除外
)

## ---- クラスタリングの既定パラメータ ----
## 03（全細胞）・07（好中球）で共通の既定値。各スクリプト内でも上書き可能。
## Elbow plot を見て変更してよいが、既定は n_dims=15 / resolution=0.6。
default_n_dims     <- 15    # PCA の使用次元数
default_resolution <- 0.6   # FindClusters の resolution

## ---- 好中球マーカー（このデータが好中球であることの確認用）----
## 生物種・対象が異なる場合は適宜変更。
neutrophil_markers <- c("S100a8", "S100a9", "Ly6g", "Itgam", "Csf3r", "Cxcr2")

## ---- ★★ 04 実行後に必ず設定 ★★ 好中球クラスタ番号 ----
## 全細胞クラスタリング(03)の結果を 04 の免疫細胞マーカー DotPlot で確認し、
## 好中球マーカー(S100a8/S100a9/Ly6g/Csf3r/Cxcr2 高発現)のクラスタ番号を列挙する。
##
## 【重要】クラスタ番号は resolution・n_dims・データにより変わるため、
##         必ず各自のデータで判定し直すこと。下記はあくまで記入例（プレースホルダ）。
##         判定できたら NULL を書き換え、05 以降を実行する。
neutrophil_clusters <- as.character(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 14))
## Non-neutrophil clusters in the final first-round annotation:
## 10 = Monocyte, 11 = NK, 13 = pDC, 15 = Basophil.

## ---- 好中球サブセット・状態を特徴づける遺伝子（08 の DotPlot 用。必要に応じて変更）----
subset_markers <- c(
  # 未成熟 -> 成熟
  "Ngp", "Ltf", "Camp", "Mmp8", "Mmp9", "Retnlg", "Cxcr4", "Sell", "Icam1",
  # 抑制性・腫瘍随伴（PMN-MDSC様）
  "Arg1", "Cd84", "Wfdc17", "Trem1",
  # I型インターフェロン応答
  "Isg15", "Irf7", "Ifit1", "Ifit3",
  # 増殖
  "Mki67", "Top2a"
)

## ---- 出力先 ----
dir_rds    <- file.path(proj_dir, "rds")
dir_output <- file.path(proj_dir, "Outputs")
dirs_needed <- c(
  dir_rds,
  file.path(dir_output, "01_QC"),
  file.path(dir_output, "02_clustering"),
  file.path(dir_output, "03_markers"),
  file.path(dir_output, "04_neutrophil_annotation"),
  file.path(dir_output, "05_neutrophil_subset_clustering"),
  file.path(dir_output, "06_neutrophil_markers"),
  file.path(dir_output, "07_neutrophil_composition"),
  file.path(dir_output, "08_neutrophil_DEG"),
  file.path(dir_output, "09_neutrophil_trajectory")
)
invisible(lapply(dirs_needed, dir.create, recursive = TRUE, showWarnings = FALSE))

## ---- 図保存ヘルパー（出力は PDF に統一）----
## filename に .png 等を渡しても .pdf に置換して保存する。
save_plot <- function(plot, filename, outdir, width = 7, height = 6) {
  filename <- sub("\\.(png|jpe?g|tiff?)$", ".pdf", filename, ignore.case = TRUE)
  path <- file.path(outdir, filename)
  ggsave(path, plot = plot, width = width, height = height, device = "pdf", bg = "white")
  message("saved: ", path)
}

message("00_setup.R 完了。sample_info を確認してください:")
print(sample_info)
