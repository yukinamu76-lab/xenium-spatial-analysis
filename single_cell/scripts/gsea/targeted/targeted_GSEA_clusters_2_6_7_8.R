#20260724 クラスター2,6,7,8に対して、全クラスターヒートマップに含まれる下記5セットでGSEAを行う




# ============================================================
# Targeted preranked GSEA
# MSigDB mouse version: 2026.1.Mm
#
# Clusters: 2, 6, 7, 8
#
# Gene sets:
#   "E2F targets",
#"Mitotic cytokinesis",
#"Hallmark interferon alpha response",
#"Interferon gamma response",
#"Regulation of granulocyte chemotaxis"
#
# Input CSVに必要な列:
#   gene
#   avg_log2FC
#
# 出力:
#   GSEA全結果CSV
#   NES行列
#   padj行列
#   global FDR行列
#   ヒートマップPNG/PDF
# ============================================================


# ------------------------------------------------------------
# 0. パッケージ
# ------------------------------------------------------------

cran_packages <- c(
  "msigdbr",
  "dplyr",
  "tidyr",
  "purrr",
  "readr",
  "ggplot2",
  "stringr",
  "tibble"
)

missing_cran <- cran_packages[
  !cran_packages %in% rownames(installed.packages())
]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install("fgsea")
}


library(msigdbr)
library(fgsea)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)
library(stringr)
library(tibble)


# ------------------------------------------------------------
# 1. 入力フォルダ
# ------------------------------------------------------------

input_dir <- file.path(getwd(), "Outputs", "08_neutrophil_DEG")


# ------------------------------------------------------------
# 2. 各クラスターのランキングCSV
# ------------------------------------------------------------

# 実際のファイル名に合わせて、ここを変更してください
input_files <- c(
  "2" = file.path(
    input_dir,
    "cluster2_full_gene_ranking.csv"
  ),
  "6" = file.path(
    input_dir,
    "cluster6_full_gene_ranking.csv"
  ),
  "7" = file.path(
    input_dir,
    "cluster7_full_gene_ranking.csv"
  ),
  "8" = file.path(
    input_dir,
    "cluster8_full_gene_ranking.csv"
  )
)


# 入力ファイル確認
file_check <- tibble(
  cluster = names(input_files),
  file = unname(input_files),
  exists = file.exists(input_files)
)

print(file_check)

if (!all(file_check$exists)) {
  
  missing_files <- file_check %>%
    filter(!exists) %>%
    pull(file)
  
  stop(
    paste0(
      "以下のファイルが見つかりません。\n",
      "ファイル名またはinput_dirを確認してください。\n\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


# ------------------------------------------------------------
# 3. 出力フォルダ
# ------------------------------------------------------------

outdir <- file.path(getwd(), "Outputs", "11_targeted_GSEA_5sets_4clusters")

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat(
  "\nOutput directory:\n",
  outdir,
  "\n\n"
)


# ------------------------------------------------------------
# 4. 使用可能なMSigDBコレクションを保存
# ------------------------------------------------------------

collections_mm <- msigdbr_collections(
  db_species = "MM"
)

print(collections_mm)

write_csv(
  collections_mm,
  file.path(
    outdir,
    "MSigDB_mouse_collections.csv"
  )
)


# ------------------------------------------------------------
# 5. MSigDB遺伝子セットの取得
# ------------------------------------------------------------

# MSigDB 2026.1.Mmでは、
# Hallmark = MH
# Reactome = M2 / CP:REACTOME
# GO BP = M5 / GO:BP


# Hallmark
hallmark_mm <- msigdbr(
  db_species = "MM",
  species = "Mus musculus",
  collection = "MH"
)


# Reactome
reactome_mm <- msigdbr(
  db_species = "MM",
  species = "Mus musculus",
  collection = "M2",
  subcollection = "CP:REACTOME"
)


# GO Biological Process
go_bp_mm <- msigdbr(
  db_species = "MM",
  species = "Mus musculus",
  collection = "M5",
  subcollection = "GO:BP"
)


cat(
  "Hallmark rows:",
  nrow(hallmark_mm),
  "\n"
)

cat(
  "Reactome rows:",
  nrow(reactome_mm),
  "\n"
)

cat(
  "GO BP rows:",
  nrow(go_bp_mm),
  "\n\n"
)


# ------------------------------------------------------------
# 6. 対象となる5遺伝子セット
# ------------------------------------------------------------

target_names <- c(
  "HALLMARK_E2F_TARGETS",
  "GOBP_MITOTIC_CYTOKINESIS",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "GOBP_REGULATION_OF_GRANULOCYTE_CHEMOTAXIS"
)


# ------------------------------------------------------------
# 7. 実際の遺伝子セット名を確認
# ------------------------------------------------------------

all_selected_collections <- bind_rows(
  hallmark_mm,
  reactome_mm,
  go_bp_mm
)


# 指定した遺伝子セットを抽出
target_msig <- all_selected_collections %>%
  filter(
    gs_name %in% target_names
  ) %>%
  distinct(
    gs_name,
    gene_symbol,
    .keep_all = TRUE
  )


found_names <- unique(
  target_msig$gs_name
)

missing_names <- setdiff(
  target_names,
  found_names
)

cat("取得できたgene set:\n")
print(found_names)


# 指定名が見つからない場合は候補名を表示して停止
if (length(missing_names) > 0) {
  
  cat("\n取得できなかったgene set:\n")
  print(missing_names)
  
  cat(
    "\nE2Fを含む候補:\n"
  )
  
  print(
    hallmark_mm %>%
      filter(
        str_detect(
          gs_name,
          regex(
            "E2F",
            ignore_case = TRUE
          )
        )
      ) %>%
      distinct(gs_name) %>%
      arrange(gs_name)
  )
  
  cat(
    "\nINTERFERONを含む候補:\n"
  )
  
  print(
    hallmark_mm %>%
      filter(
        str_detect(
          gs_name,
          regex(
            "INTERFERON",
            ignore_case = TRUE
          )
        )
      ) %>%
      distinct(gs_name) %>%
      arrange(gs_name)
  )
  
  cat(
    "\nMITOTICまたはSPINDLEを含むReactome候補:\n"
  )
  
  print(
    reactome_mm %>%
      filter(
        str_detect(
          gs_name,
          regex(
            "MITOTIC|SPINDLE",
            ignore_case = TRUE
          )
        )
      ) %>%
      distinct(gs_name) %>%
      arrange(gs_name)
  )
  
  cat(
    "\nCELL_CHEMOTAXISを含むGO BP候補:\n"
  )
  
  print(
    go_bp_mm %>%
      filter(
        str_detect(
          gs_name,
          regex(
            "CELL.*CHEMOTAXIS|CHEMOTAXIS",
            ignore_case = TRUE
          )
        )
      ) %>%
      distinct(gs_name) %>%
      arrange(gs_name)
  )
  
  stop(
    paste0(
      "指定した遺伝子セットがすべて見つかりませんでした: ",
      paste(
        missing_names,
        collapse = ", "
      )
    )
  )
}


# 5セットすべて取得できたことを確認
if (length(found_names) != length(target_names)) {
  stop(
    "対象gene setの数が5個ではありません。"
  )
}


# ------------------------------------------------------------
# 8. Gene set情報を保存
# ------------------------------------------------------------

gene_set_summary <- target_msig %>%
  count(
    gs_name,
    name = "n_genes_in_MSigDB"
  ) %>%
  mutate(
    gene_set_order = match(
      gs_name,
      target_names
    )
  ) %>%
  arrange(gene_set_order) %>%
  select(
    -gene_set_order
  )

print(gene_set_summary)

write_csv(
  gene_set_summary,
  file.path(
    outdir,
    "selected_gene_set_summary.csv"
  )
)


write_csv(
  target_msig %>%
    select(
      gene_set = gs_name,
      gene = gene_symbol,
      everything()
    ),
  file.path(
    outdir,
    "selected_gene_set_members.csv"
  )
)


# ------------------------------------------------------------
# 9. fgsea用pathwayリスト
# ------------------------------------------------------------

pathways <- split(
  target_msig$gene_symbol,
  target_msig$gs_name
)

pathways <- lapply(
  pathways,
  unique
)


# target_namesの順番に並べる
pathways <- pathways[
  target_names
]


cat(
  "\n各gene setの遺伝子数:\n"
)

print(
  lengths(pathways)
)


# ------------------------------------------------------------
# 10. ランキングデータ作成関数
# ------------------------------------------------------------

prepare_ranks <- function(
    file_path,
    cluster_id
) {
  
  message(
    "\nReading cluster ",
    cluster_id,
    ": ",
    file_path
  )
  
  dat <- read_csv(
    file_path,
    show_col_types = FALSE
  )
  
  
  cat(
    "Cluster ",
    cluster_id,
    " columns:\n",
    sep = ""
  )
  
  print(
    colnames(dat)
  )
  
  
  required_columns <- c(
    "gene",
    "avg_log2FC"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(dat)
  )
  
  
  if (length(missing_columns) > 0) {
    
    stop(
      paste0(
        "Cluster ",
        cluster_id,
        "のCSVに必要な列がありません: ",
        paste(
          missing_columns,
          collapse = ", "
        ),
        "\n現在の列名: ",
        paste(
          colnames(dat),
          collapse = ", "
        )
      )
    )
  }
  
  
  rank_df <- dat %>%
    transmute(
      gene = as.character(gene),
      avg_log2FC = suppressWarnings(
        as.numeric(avg_log2FC)
      )
    ) %>%
    filter(
      !is.na(gene),
      gene != "",
      !is.na(avg_log2FC),
      is.finite(avg_log2FC)
    ) %>%
    
    # 同じ遺伝子が複数行ある場合、
    # |avg_log2FC|が最大の行を使用
    group_by(gene) %>%
    slice_max(
      order_by = abs(avg_log2FC),
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup() %>%
    arrange(
      desc(avg_log2FC)
    )
  
  
  ranks <- rank_df$avg_log2FC
  names(ranks) <- rank_df$gene
  
  ranks <- sort(
    ranks,
    decreasing = TRUE
  )
  
  
  if (anyDuplicated(names(ranks)) > 0) {
    stop(
      paste0(
        "Cluster ",
        cluster_id,
        "のランキングに重複遺伝子があります。"
      )
    )
  }
  
  
  if (length(ranks) < 1000) {
    warning(
      paste0(
        "Cluster ",
        cluster_id,
        "のランキング遺伝子数が少ないです: ",
        length(ranks),
        " genes"
      )
    )
  }
  
  
  cat(
    "Cluster ",
    cluster_id,
    " ranking genes: ",
    length(ranks),
    "\n",
    sep = ""
  )
  
  cat(
    "Positive statistics: ",
    sum(ranks > 0),
    "\n",
    sep = ""
  )
  
  cat(
    "Negative statistics: ",
    sum(ranks < 0),
    "\n",
    sep = ""
  )
  
  cat(
    "Zero statistics: ",
    sum(ranks == 0),
    "\n",
    sep = ""
  )
  
  
  list(
    ranks = ranks,
    rank_df = rank_df
  )
}


# ------------------------------------------------------------
# 11. fgsea実行関数
# ------------------------------------------------------------

run_targeted_fgsea <- function(
    file_path,
    cluster_id,
    pathways,
    outdir
) {
  
  prepared <- prepare_ranks(
    file_path = file_path,
    cluster_id = cluster_id
  )
  
  ranks <- prepared$ranks
  rank_df <- prepared$rank_df
  
  
  # 各gene setとランキングとの重複数
  overlap_df <- tibble(
    pathway = names(pathways),
    n_genes_original = lengths(pathways),
    n_genes_in_ranking = map_int(
      pathways,
      ~ sum(.x %in% names(ranks))
    )
  )
  
  
  print(
    overlap_df
  )
  
  
  write_csv(
    overlap_df,
    file.path(
      outdir,
      paste0(
        "cluster",
        cluster_id,
        "_gene_set_overlap.csv"
      )
    )
  )
  
  
  # 実際に使用したランキングを保存
  write_csv(
    rank_df,
    file.path(
      outdir,
      paste0(
        "cluster",
        cluster_id,
        "_ranking_used.csv"
      )
    )
  )
  
  
  set.seed(1234)
  
  
  fgsea_res <- fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 0,
    nPermSimple = 100000
  )
  
  
  fgsea_res <- as.data.frame(
    fgsea_res
  ) %>%
    as_tibble()
  
  
  # leadingEdgeを文字列へ変換
  if ("leadingEdge" %in% colnames(fgsea_res)) {
    
    fgsea_res <- fgsea_res %>%
      mutate(
        leadingEdge = map_chr(
          leadingEdge,
          ~ paste(
            .x,
            collapse = ";"
          )
        )
      )
    
  } else {
    
    fgsea_res <- fgsea_res %>%
      mutate(
        leadingEdge = NA_character_
      )
  }
  
  
  fgsea_res <- fgsea_res %>%
    mutate(
      cluster = as.character(cluster_id)
    ) %>%
    select(
      cluster,
      pathway,
      size,
      ES,
      NES,
      pval,
      padj,
      log2err,
      leadingEdge
    )
  
  
  # clusterごとの結果も保存
  write_csv(
    fgsea_res,
    file.path(
      outdir,
      paste0(
        "cluster",
        cluster_id,
        "_targeted_GSEA_results.csv"
      )
    )
  )
  
  
  fgsea_res
}


# ------------------------------------------------------------
# 12. 4クラスターでGSEAを実行
# ------------------------------------------------------------

gsea_all <- imap_dfr(
  input_files,
  ~ run_targeted_fgsea(
    file_path = .x,
    cluster_id = .y,
    pathways = pathways,
    outdir = outdir
  )
)


# ------------------------------------------------------------
# 13. 5セット × 4クラスターのglobal FDR
# ------------------------------------------------------------

# 各クラスター内のpadj:
# そのクラスターの5セットを対象にBH補正
#
# padj_global_20tests:
# 5セット × 4クラスター = 20検定をまとめてBH補正

gsea_all <- gsea_all %>%
  mutate(
    padj_global_20tests = p.adjust(
      pval,
      method = "BH"
    )
  )


# ------------------------------------------------------------
# 14. 表示名
# ------------------------------------------------------------

pathway_labels <- c(
  "HALLMARK_E2F_TARGETS" =
    "E2F targets",
  
  "GOBP_MITOTIC_CYTOKINESIS" =
    "Mitotic cytokinesis",
  
  "HALLMARK_INTERFERON_ALPHA_RESPONSE" =
    "Hallmark interferon alpha response",
  
  "HALLMARK_INTERFERON_GAMMA_RESPONSE" =
    "Interferon gamma response",
  
  "GOBP_REGULATION_OF_GRANULOCYTE_CHEMOTAXIS" =
    "Regulation of granulocyte chemotaxis"
)


gsea_all <- gsea_all %>%
  mutate(
    pathway_label = unname(
      pathway_labels[pathway]
    )
  )


if (any(is.na(gsea_all$pathway_label))) {
  warning(
    "表示名に変換できなかったpathwayがあります。"
  )
}


# ------------------------------------------------------------
# 15. 全結果CSV保存
# ------------------------------------------------------------

write_csv(
  gsea_all,
  file.path(
    outdir,
    "targeted_GSEA_5sets_4clusters_full_results.csv"
  )
)


heatmap_input <- gsea_all %>%
  select(
    pathway,
    pathway_label,
    cluster,
    size,
    ES,
    NES,
    pval,
    padj,
    padj_global_20tests,
    leadingEdge
  )


write_csv(
  heatmap_input,
  file.path(
    outdir,
    "targeted_GSEA_5sets_4clusters_heatmap_input.csv"
  )
)


# ------------------------------------------------------------
# 16. NES行列
# ------------------------------------------------------------

nes_wide <- heatmap_input %>%
  select(
    pathway_label,
    cluster,
    NES
  ) %>%
  pivot_wider(
    names_from = cluster,
    values_from = NES,
    names_prefix = "cluster"
  )


write_csv(
  nes_wide,
  file.path(
    outdir,
    "targeted_GSEA_5sets_4clusters_NES_matrix.csv"
  )
)


# ------------------------------------------------------------
# 17. padj行列
# ------------------------------------------------------------

padj_wide <- heatmap_input %>%
  select(
    pathway_label,
    cluster,
    padj
  ) %>%
  pivot_wider(
    names_from = cluster,
    values_from = padj,
    names_prefix = "cluster"
  )


write_csv(
  padj_wide,
  file.path(
    outdir,
    "targeted_GSEA_5sets_4clusters_padj_matrix.csv"
  )
)


# ------------------------------------------------------------
# 18. global FDR行列
# ------------------------------------------------------------

global_padj_wide <- heatmap_input %>%
  select(
    pathway_label,
    cluster,
    padj_global_20tests
  ) %>%
  pivot_wider(
    names_from = cluster,
    values_from = padj_global_20tests,
    names_prefix = "cluster"
  )


write_csv(
  global_padj_wide,
  file.path(
    outdir,
    "targeted_GSEA_5sets_4clusters_global_padj_matrix.csv"
  )
)


# ------------------------------------------------------------
# 19. ヒートマップの順番
# ------------------------------------------------------------

# Slingshotのpseudotime順
cluster_order <- c(
  "8",
  "6",
  "7",
  "2"
)


pathway_order <- c(
  "E2F targets",
  "Mitotic cytokinesis",
  "Hallmark interferon alpha response",
  "Interferon gamma response",
  "Regulation of granulocyte chemotaxis"
)


plot_df <- gsea_all %>%
  mutate(
    cluster = factor(
      cluster,
      levels = cluster_order
    ),
    
    # ggplotでは下から上に描画されるためrev()
    pathway_label = factor(
      pathway_label,
      levels = rev(pathway_order)
    )
  )


# ------------------------------------------------------------
# 20. セル内の表示
# ------------------------------------------------------------

plot_df <- plot_df %>%
  mutate(
    label = case_when(
      
      is.na(NES) ~
        "NA\n(padj=NA)",
      
      is.na(padj) ~
        paste0(
          sprintf("%.2f", NES),
          "\n(padj=NA)"
        ),
      
      TRUE ~
        paste0(
          sprintf("%.2f", NES),
          "\n(padj=",
          formatC(
            padj,
            format = "e",
            digits = 1
          ),
          ")"
        )
    )
  )


# 色の範囲を正負対称にする
if (all(is.na(plot_df$NES))) {
  
  nes_limit <- 2
  
} else {
  
  nes_limit <- max(
    abs(plot_df$NES),
    na.rm = TRUE
  )
  
  nes_limit <- ceiling(
    nes_limit * 10
  ) / 10
}


# ------------------------------------------------------------
# 21. クラスター内padjを表示するヒートマップ
# ------------------------------------------------------------

p_heatmap <- ggplot(
  plot_df,
  aes(
    x = cluster,
    y = pathway_label,
    fill = NES
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -nes_limit,
      nes_limit
    ),
    na.value = "grey65",
    name = "NES"
  ) +
  
  geom_text(
    aes(
      label = label
    ),
    size = 3.8,
    lineheight = 0.95,
    color = "black"
  ) +
  
  scale_x_discrete(
    drop = FALSE
  ) +
  
  scale_y_discrete(
    drop = FALSE
  ) +
  
  labs(
    x = "Cluster",
    y = NULL,
    title = NULL
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    axis.text.x = element_text(
      size = 13,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 12,
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 14,
      color = "black"
    ),
    
    axis.ticks = element_blank(),
    
    panel.grid = element_blank(),
    
    legend.title = element_text(
      size = 13
    ),
    
    legend.text = element_text(
      size = 11
    ),
    
    plot.margin = margin(
      t = 10,
      r = 20,
      b = 10,
      l = 10
    )
  )


print(
  p_heatmap
)


# ------------------------------------------------------------
# 22. ヒートマップ保存
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    outdir,
    "targeted_GSEA_5sets_4clusters_heatmap.png"
  ),
  plot = p_heatmap,
  width = 10.5,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# cairo_pdfが使えない環境にも対応
tryCatch(
  
  ggsave(
    filename = file.path(
      outdir,
      "targeted_GSEA_5sets_4clusters_heatmap.pdf"
    ),
    plot = p_heatmap,
    width = 10.5,
    height = 6.5,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  ),
  
  error = function(e) {
    
    message(
      "cairo_pdfで保存できなかったため、標準PDFを使用します。"
    )
    
    ggsave(
      filename = file.path(
        outdir,
        "targeted_GSEA_5sets_4clusters_heatmap.pdf"
      ),
      plot = p_heatmap,
      width = 10.5,
      height = 6.5,
      units = "in",
      device = "pdf",
      bg = "white"
    )
  }
)


# ------------------------------------------------------------
# 23. global FDRを表示するヒートマップ
# ------------------------------------------------------------

plot_df_global <- gsea_all %>%
  mutate(
    cluster = factor(
      cluster,
      levels = cluster_order
    ),
    
    pathway_label = factor(
      pathway_label,
      levels = rev(pathway_order)
    ),
    
    label = case_when(
      
      is.na(NES) ~
        "NA\n(FDR=NA)",
      
      is.na(padj_global_20tests) ~
        paste0(
          sprintf("%.2f", NES),
          "\n(FDR=NA)"
        ),
      
      TRUE ~
        paste0(
          sprintf("%.2f", NES),
          "\n(FDR=",
          formatC(
            padj_global_20tests,
            format = "e",
            digits = 1
          ),
          ")"
        )
    )
  )


p_heatmap_global <- ggplot(
  plot_df_global,
  aes(
    x = cluster,
    y = pathway_label,
    fill = NES
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -nes_limit,
      nes_limit
    ),
    na.value = "grey65",
    name = "NES"
  ) +
  
  geom_text(
    aes(
      label = label
    ),
    size = 3.8,
    lineheight = 0.95,
    color = "black"
  ) +
  
  scale_x_discrete(
    drop = FALSE
  ) +
  
  scale_y_discrete(
    drop = FALSE
  ) +
  
  labs(
    x = "Cluster",
    y = NULL,
    title = NULL
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    axis.text.x = element_text(
      size = 13,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 12,
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 14,
      color = "black"
    ),
    
    axis.ticks = element_blank(),
    
    panel.grid = element_blank(),
    
    legend.title = element_text(
      size = 13
    ),
    
    legend.text = element_text(
      size = 11
    ),
    
    plot.margin = margin(
      t = 10,
      r = 20,
      b = 10,
      l = 10
    )
  )


print(
  p_heatmap_global
)


ggsave(
  filename = file.path(
    outdir,
    "targeted_GSEA_5sets_4clusters_heatmap_globalFDR.png"
  ),
  plot = p_heatmap_global,
  width = 10.5,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


tryCatch(
  
  ggsave(
    filename = file.path(
      outdir,
      "targeted_GSEA_5sets_4clusters_heatmap_globalFDR.pdf"
    ),
    plot = p_heatmap_global,
    width = 10.5,
    height = 6.5,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  ),
  
  error = function(e) {
    
    message(
      "cairo_pdfで保存できなかったため、標準PDFを使用します。"
    )
    
    ggsave(
      filename = file.path(
        outdir,
        "targeted_GSEA_5sets_4clusters_heatmap_globalFDR.pdf"
      ),
      plot = p_heatmap_global,
      width = 10.5,
      height = 6.5,
      units = "in",
      device = "pdf",
      bg = "white"
    )
  }
)


# ------------------------------------------------------------
# 24. 最終結果をコンソールに表示
# ------------------------------------------------------------

final_table <- gsea_all %>%
  mutate(
    cluster_sort = factor(
      cluster,
      levels = cluster_order
    ),
    
    pathway_sort = factor(
      pathway_label,
      levels = pathway_order
    )
  ) %>%
  arrange(
    pathway_sort,
    cluster_sort
  ) %>%
  select(
    cluster,
    pathway,
    pathway_label,
    size,
    ES,
    NES,
    pval,
    padj,
    padj_global_20tests,
    leadingEdge
  )


print(
  final_table,
  n = Inf,
  width = Inf
)


# ------------------------------------------------------------
# Global FDR一覧を保存
# ------------------------------------------------------------

global_fdr_table <- gsea_all %>%
  select(
    cluster,
    pathway,
    pathway_label,
    NES,
    pval,
    padj,
    padj_global_20tests
  ) %>%
  arrange(
    pathway_label,
    factor(cluster, levels = c("8","6","7","2"))
  )

write_csv(
  global_fdr_table,
  file.path(
    outdir,
    "targeted_GSEA_5sets_4clusters_globalFDR_results.csv"
  )
)


# ------------------------------------------------------------
# Global FDR matrix
# ------------------------------------------------------------

global_fdr_matrix <- gsea_all %>%
  select(
    pathway_label,
    cluster,
    padj_global_20tests
  ) %>%
  pivot_wider(
    names_from = cluster,
    values_from = padj_global_20tests,
    names_prefix = "cluster"
  )

write_csv(
  global_fdr_matrix,
  file.path(
    outdir,
    "targeted_GSEA_5sets_4clusters_globalFDR_matrix.csv"
  )
)

cat(
  "\n========================================\n",
  "Targeted GSEA completed.\n",
  "MSigDB: 2026.1.Mm\n",
  "Output directory:\n",
  outdir,
  "\n========================================\n"
)
