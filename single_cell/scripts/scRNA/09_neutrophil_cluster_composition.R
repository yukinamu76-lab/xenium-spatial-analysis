## ============================================================
## 09_neutrophil_cluster_composition.R
## Manuscript outputs: sample-level cluster composition for Table EV6 and source data for Fig. EV4D
##   → 好中球クラスタごとの群間 組成比較（③ 2回目UMAPの後に実行）
##
##   「どのクラスタが群間で増える/減るか」を、サンプル単位のクラスタ比率で比較する。
##   出力は csv（比率表・検定結果）と pdf（積み上げ棒・箱ひげ・log2FC バー）。
##
## 【注意】サンプル数が少ない場合（例: 各群2）、Wilcoxon/t 検定はあくまで記述的な
##         参考値。可能なら speckle::propeller（レプリケート考慮）を併用する。
## ============================================================
source(file.path("scripts", "scRNA", "00_setup.R"))

comp_dir <- file.path(dir_output, "07_neutrophil_composition")
neut     <- readRDS(file.path(dir_rds, "07_neutrophil_clustered.rds"))

## 群ラベル（00_setup.R の group_levels より。先頭=基準群）
g_ref  <- group_levels[1]   # 基準群（対照）
g_test <- group_levels[2]   # 比較群

meta <- neut@meta.data %>%
  as_tibble(rownames = "cell") %>%
  mutate(cluster = seurat_clusters)

## ---- サンプルごとのクラスタ比率 ----
prop_df <- meta %>%
  count(sample_id, group, cluster, name = "n_cells") %>%
  group_by(sample_id) %>%
  mutate(total_cells = sum(n_cells), prop = n_cells / total_cells) %>%
  ungroup()
write_csv(prop_df, file.path(comp_dir, "cluster_proportions_per_sample.csv"))

## ---- 積み上げ棒グラフ: サンプルごとのクラスタ構成 ----
p_stack <- ggplot(prop_df, aes(x = sample_id, y = prop, fill = cluster)) +
  geom_col() +
  facet_grid(~group, scales = "free_x", space = "free_x") +
  labs(y = "Proportion of cells", x = NULL, fill = "Cluster") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_stack, "01_stacked_bar_per_sample.pdf", comp_dir, width = 8, height = 6)

## ---- クラスタごとの比率比較（群平均±サンプル点）----
p_box <- ggplot(prop_df, aes(x = group, y = prop, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2) +
  facet_wrap(~cluster, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  labs(y = "Proportion of cells", x = NULL) +
  theme_minimal(base_size = 11)
save_plot(p_box, "02_proportion_boxplot_per_cluster.pdf", comp_dir, width = 10, height = 8)

## ---- 統計検定: クラスタごとに 比較群 vs 基準群 の比率差を検定 ----
cluster_test <- prop_df %>%
  group_by(cluster) %>%
  summarise(
    mean_ref  = mean(prop[group == g_ref]),
    mean_test = mean(prop[group == g_test]),
    log2FC_test_vs_ref = log2((mean_test + 1e-4) / (mean_ref + 1e-4)),
    p_wilcox = tryCatch(
      wilcox.test(prop[group == g_test], prop[group == g_ref])$p.value,
      error = function(e) NA_real_),
    p_ttest = tryCatch(
      t.test(prop[group == g_test], prop[group == g_ref])$p.value,
      error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  rename(!!paste0("mean_", g_ref)  := mean_ref,
         !!paste0("mean_", g_test) := mean_test) %>%
  arrange(desc(abs(log2FC_test_vs_ref)))
write_csv(cluster_test, file.path(comp_dir, "cluster_composition_test.csv"))
print(cluster_test)

## ---- speckle::propeller（利用可能ならレプリケートを考慮した検定として併用）----
if (requireNamespace("speckle", quietly = TRUE)) {
  library(speckle)
  prop_test <- tryCatch(
    propeller(clusters = meta$cluster, sample = meta$sample_id, group = meta$group),
    error = function(e) { message("propeller 実行エラー: ", e$message); NULL }
  )
  if (!is.null(prop_test)) {
    write_csv(as_tibble(prop_test, rownames = "cluster"),
              file.path(comp_dir, "cluster_composition_propeller.csv"))
    print(prop_test)
  }
} else {
  message("speckle が無いため propeller はスキップ（Wilcoxon/t 検定のみ）")
}

## ---- 増加/減少クラスタの可視化 (log2FC バープロット) ----
lab_up   <- sprintf("%s で増加", g_test)
lab_down <- sprintf("%s で増加", g_ref)
p_fc <- ggplot(cluster_test, aes(x = reorder(cluster, log2FC_test_vs_ref),
                                 y = log2FC_test_vs_ref,
                                 fill = log2FC_test_vs_ref > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
                    labels = c("TRUE" = lab_up, "FALSE" = lab_down), name = NULL) +
  labs(x = "Cluster",
       y = sprintf("log2FC (%s vs %s, クラスタ比率)", g_test, g_ref)) +
  theme_minimal(base_size = 12)
save_plot(p_fc, "03_cluster_log2FC_barplot.pdf", comp_dir, width = 7, height = 6)

message("09_neutrophil_cluster_composition.R 完了")
message(sprintf("%s/07_neutrophil_composition/cluster_composition_test.csv で増減クラスタを確認してください", dir_output))
