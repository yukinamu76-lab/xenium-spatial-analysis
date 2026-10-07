# Manuscript outputs: comprehensive cluster 5 GSEA contributing to Dataset EV7 and Fig. EV4E
# Comprehensive pre-ranked GSEA for neutrophil cluster 5 vs retained clusters
# クラスター5vs他クラスターのlog2FC順の全遺伝子に対してGSEAを行う
#  インプットは、cluster5_full_gene_ranking.csv

# ============================================================
# Cluster 5 vs rest: preranked GSEA
# Collections:
#   1. MSigDB Hallmark
#   2. GO Biological Process
#   3. Reactome
#
# Ranking:
#   avg_log2FC
#   positive NES = cluster 5で高い
#   negative NES = cluster 5で低い（restで高い）
# ============================================================


# ------------------------------------------------------------
# 0. 必要パッケージ
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install("fgsea")
}

cran_packages <- c(
  "msigdbr",
  "dplyr",
  "ggplot2",
  "stringr",
  "purrr",
  "tibble"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
library(tibble)


# ------------------------------------------------------------
# 1. 入出力パス
# ------------------------------------------------------------

input_file <- file.path(getwd(), "Outputs", "08_neutrophil_DEG", "cluster5_full_gene_ranking.csv")

outdir <- file.path(getwd(), "Outputs", "10_GSEA", "cluster5")

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. DEGファイル読み込み
# ------------------------------------------------------------

deg <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 列名確認
print(colnames(deg))
print(dim(deg))
print(head(deg))

required_cols <- c(
  "gene",
  "avg_log2FC"
)

missing_cols <- setdiff(
  required_cols,
  colnames(deg)
)

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "必要な列がありません: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 3. fgsea用ランキング作成
# ------------------------------------------------------------

rank_df <- deg %>%
  transmute(
    gene = as.character(gene),
    avg_log2FC = as.numeric(avg_log2FC)
  ) %>%
  filter(
    !is.na(gene),
    gene != "",
    !is.na(avg_log2FC),
    is.finite(avg_log2FC)
  ) %>%
  group_by(gene) %>%
  slice_max(
    order_by = abs(avg_log2FC),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(
    desc(avg_log2FC),
    gene
  )

ranks <- rank_df$avg_log2FC
names(ranks) <- rank_df$gene

ranks <- sort(
  ranks,
  decreasing = TRUE
)

cat("ランキング遺伝子数:", length(ranks), "\n")
cat("重複遺伝子数:", anyDuplicated(names(ranks)), "\n")
cat("正のlog2FC:", sum(ranks > 0), "\n")
cat("負のlog2FC:", sum(ranks < 0), "\n")

print(head(ranks, 10))
print(tail(ranks, 10))


# ランキング入力値を保存
rank_input <- data.frame(
  rank = seq_along(ranks),
  gene = names(ranks),
  avg_log2FC = as.numeric(ranks),
  direction = ifelse(
    ranks > 0,
    "Higher_in_cluster5",
    ifelse(
      ranks < 0,
      "Lower_in_cluster5",
      "No_change"
    )
  ),
  check.names = FALSE
)

write.csv(
  rank_input,
  file.path(
    outdir,
    "cluster5_GSEA_ranked_gene_input.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. msigdbrコレクション確認
# ------------------------------------------------------------

# 現在利用可能なマウスコレクションを表示
collection_info <- tryCatch(
  msigdbr_collections(
    db_species = "MM"
  ),
  error = function(e) {
    msigdbr_collections()
  }
)

print(collection_info)

write.csv(
  collection_info,
  file.path(
    outdir,
    "msigdbr_available_collections.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. MSigDB遺伝子セット取得用関数
#    新旧msigdbr両方に対応
# ------------------------------------------------------------

get_msig_mouse <- function(
    collection_name,
    subcollection_name = NULL
) {
  
  # 新しい書式を優先
  result <- tryCatch(
    {
      msigdbr(
        db_species = "MM",
        species = "Mus musculus",
        collection = collection_name,
        subcollection = subcollection_name
      )
    },
    error = function(e1) {
      
      message(
        "新しいmsigdbr書式で取得できなかったため、",
        "旧書式を試します。"
      )
      
      msigdbr(
        species = "Mus musculus",
        category = collection_name,
        subcategory = subcollection_name
      )
    }
  )
  
  return(result)
}


# ------------------------------------------------------------
# 6. Hallmark、GO BP、Reactomeを取得
# ------------------------------------------------------------

msig_hallmark <- get_msig_mouse(
  collection_name = "H"
)

msig_go_bp <- get_msig_mouse(
  collection_name = "C5",
  subcollection_name = "GO:BP"
)

msig_reactome <- get_msig_mouse(
  collection_name = "C2",
  subcollection_name = "CP:REACTOME"
)


# ------------------------------------------------------------
# 7. pathwayリストへ変換
# ------------------------------------------------------------

make_pathway_list <- function(msig_df) {
  
  if (!all(c("gs_name", "gene_symbol") %in% colnames(msig_df))) {
    stop(
      "msigdbr出力にgs_nameまたはgene_symbol列がありません。"
    )
  }
  
  pathways <- split(
    x = msig_df$gene_symbol,
    f = msig_df$gs_name
  )
  
  pathways <- lapply(
    pathways,
    function(x) unique(
      x[!is.na(x) & x != ""]
    )
  )
  
  return(pathways)
}

pathways_hallmark <- make_pathway_list(
  msig_hallmark
)

pathways_go_bp <- make_pathway_list(
  msig_go_bp
)

pathways_reactome <- make_pathway_list(
  msig_reactome
)

cat(
  "Hallmark pathways:",
  length(pathways_hallmark),
  "\n"
)

cat(
  "GO BP pathways:",
  length(pathways_go_bp),
  "\n"
)

cat(
  "Reactome pathways:",
  length(pathways_reactome),
  "\n"
)



# GSEA実行部分

# ------------------------------------------------------------
# 8. fgsea実行関数
# ------------------------------------------------------------

run_fgsea_collection <- function(
    pathways,
    ranks,
    collection_label,
    min_size = 10,
    max_size = 1000
) {
  
  # ランキング中に存在する遺伝子だけに絞ったpathwayサイズを確認
  pathway_overlap <- tibble(
    pathway = names(pathways),
    original_size = lengths(pathways),
    overlap_size = vapply(
      pathways,
      function(x) {
        sum(x %in% names(ranks))
      },
      numeric(1)
    )
  )
  
  write.csv(
    pathway_overlap,
    file.path(
      outdir,
      paste0(
        "cluster5_",
        collection_label,
        "_pathway_overlap.csv"
      )
    ),
    row.names = FALSE
  )
  
  set.seed(123)
  
  res <- fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = min_size,
    maxSize = max_size,
    eps = 0
  )
  
  res_df <- as.data.frame(res) %>%
    arrange(
      padj,
      desc(abs(NES))
    ) %>%
    mutate(
      collection = collection_label,
      enrichment_direction = case_when(
        NES > 0 ~ "Higher_in_cluster5",
        NES < 0 ~ "Lower_in_cluster5",
        TRUE ~ "Neutral"
      )
    )
  
  return(res_df)
}


# ------------------------------------------------------------
# 9. 3コレクションでGSEA実行
# ------------------------------------------------------------

fgsea_hallmark <- run_fgsea_collection(
  pathways = pathways_hallmark,
  ranks = ranks,
  collection_label = "Hallmark",
  min_size = 10,
  max_size = 500
)

fgsea_go_bp <- run_fgsea_collection(
  pathways = pathways_go_bp,
  ranks = ranks,
  collection_label = "GO_BP",
  min_size = 10,
  max_size = 1000
)

fgsea_reactome <- run_fgsea_collection(
  pathways = pathways_reactome,
  ranks = ranks,
  collection_label = "Reactome",
  min_size = 10,
  max_size = 1000
)


# 統合結果
fgsea_all <- bind_rows(
  fgsea_hallmark,
  fgsea_go_bp,
  fgsea_reactome
)


# 結果のCSV保存

# ------------------------------------------------------------
# 10. CSV保存用関数
# ------------------------------------------------------------

prepare_fgsea_for_csv <- function(res_df) {
  
  res_df %>%
    mutate(
      leadingEdge = vapply(
        leadingEdge,
        function(x) {
          paste(x, collapse = ";")
        },
        FUN.VALUE = character(1)
      )
    )
}


fgsea_hallmark_save <- prepare_fgsea_for_csv(
  fgsea_hallmark
)

fgsea_go_bp_save <- prepare_fgsea_for_csv(
  fgsea_go_bp
)

fgsea_reactome_save <- prepare_fgsea_for_csv(
  fgsea_reactome
)

fgsea_all_save <- prepare_fgsea_for_csv(
  fgsea_all
)


write.csv(
  fgsea_hallmark_save,
  file.path(
    outdir,
    "cluster5_GSEA_Hallmark_all_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  fgsea_go_bp_save,
  file.path(
    outdir,
    "cluster5_GSEA_GO_BP_all_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  fgsea_reactome_save,
  file.path(
    outdir,
    "cluster5_GSEA_Reactome_all_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  fgsea_all_save,
  file.path(
    outdir,
    "cluster5_GSEA_all_collections_results.csv"
  ),
  row.names = FALSE
)


# FDR < 0.05のみ
fgsea_significant <- fgsea_all_save %>%
  filter(
    !is.na(padj),
    padj < 0.05
  ) %>%
  arrange(
    collection,
    padj,
    desc(abs(NES))
  )

write.csv(
  fgsea_significant,
  file.path(
    outdir,
    "cluster5_GSEA_all_collections_FDR005.csv"
  ),
  row.names = FALSE
)

cat(
  "FDR < 0.05のpathway数:",
  nrow(fgsea_significant),
  "\n"
)
