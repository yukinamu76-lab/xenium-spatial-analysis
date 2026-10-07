suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr); library(tibble)})
source("config.R")
long_file <- file.path(output_dir, "01_all_clusters_GSEA_long.csv")
class_file <- file.path(output_dir, "02_classification_used_v3.csv")
if (!file.exists(long_file) || !file.exists(class_file)) stop("Run 01_build_all_cluster_tables.R first.")

gsea <- read_csv(long_file, show_col_types = FALSE) |>
  mutate(cluster = as.character(cluster), padj_safe = pmax(as.numeric(padj), 1e-300), NES = as.numeric(NES))
classification <- read_csv(class_file, show_col_types = FALSE)

rescale01 <- function(x) {
  x[!is.finite(x)] <- NA_real_; r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) == 0) return(rep(0, length(x)))
  replace((x-r[1])/diff(r), is.na(x), 0)
}

stats <- gsea |> group_by(pathway) |> summarise(
  min_FDR = suppressWarnings(min(padj_safe, na.rm=TRUE)),
  max_abs_NES = suppressWarnings(max(abs(NES), na.rm=TRUE)),
  dominant_cluster = {i <- which.max(replace(abs(NES), !is.finite(NES), -Inf)); if(length(i)==0) NA_character_ else cluster[i][1]},
  dominant_NES = {i <- which.max(replace(abs(NES), !is.finite(NES), -Inf)); if(length(i)==0) NA_real_ else NES[i][1]},
  n_significant_clusters = sum(padj_safe < fdr_cutoff & abs(NES) >= nes_cutoff, na.rm=TRUE),
  n_FDR_significant_clusters = sum(padj_safe < fdr_cutoff, na.rm=TRUE),
  top1_abs_NES = {z<-sort(abs(NES[is.finite(NES)]),decreasing=TRUE); if(length(z)) z[1] else NA_real_},
  top2_abs_NES = {z<-sort(abs(NES[is.finite(NES)]),decreasing=TRUE); if(length(z)>=2) z[2] else 0},
  rest_mean_abs_NES = {z<-sort(abs(NES[is.finite(NES)]),decreasing=TRUE); if(length(z)>2) mean(z[3:length(z)]) else 0},
  .groups="drop"
) |> mutate(
  min_FDR = if_else(is.infinite(min_FDR), NA_real_, min_FDR),
  max_abs_NES = if_else(is.infinite(max_abs_NES), NA_real_, max_abs_NES),
  specificity_top1 = top1_abs_NES-top2_abs_NES,
  specificity_top2 = rowMeans(cbind(top1_abs_NES,top2_abs_NES),na.rm=TRUE)-rest_mean_abs_NES,
  significance_raw=-log10(min_FDR), rarity_raw=1/pmax(n_FDR_significant_clusters,1)
) |> left_join(classification, by="pathway") |>
  filter(eligible_for_v3_selection, functional_block %in% block_order)

stats <- stats |> group_by(functional_block) |> mutate(
  significance_score=rescale01(significance_raw), effect_score=rescale01(max_abs_NES),
  specificity_score=rescale01(pmax(specificity_top1,specificity_top2,na.rm=TRUE)), rarity_score=rescale01(rarity_raw),
  interpretability_score = case_when(
    collection == "Hallmark" ~ 1,
    collection == "Reactome" ~ 0.8,
    collection == "GO BP" ~ 0.6,
    TRUE ~ 0.4
  ),
  composite_score=selection_weights[["significance"]]*significance_score +
    selection_weights[["effect"]]*effect_score + selection_weights[["specificity"]]*specificity_score +
    selection_weights[["rarity"]]*rarity_score + selection_weights[["interpretability"]]*interpretability_score
) |> ungroup()

quota_for <- function(block) {
  if(length(block)==0 || is.na(block) || !nzchar(block)) return(default_pathways_per_block)
  if(block %in% names(pathways_per_block)) return(as.integer(pathways_per_block[[block]]))
  default_pathways_per_block
}

select_block <- function(df, block_name) {
  block <- as.character(block_name)[1]; quota <- min(quota_for(block), nrow(df))
  if(quota<=0 || nrow(df)==0) return(df[0,,drop=FALSE])
  eligible <- df |> filter(is.finite(composite_score), min_FDR < fdr_cutoff, max_abs_NES >= nes_cutoff) |>
    arrange(desc(composite_score), min_FDR, desc(max_abs_NES))
  if(nrow(eligible)==0) eligible <- df |> filter(is.finite(composite_score)) |> arrange(desc(composite_score))
  first <- eligible |> mutate(dc=replace_na(dominant_cluster,"NA")) |> group_by(dc) |>
    slice_head(n=max_per_dominant_cluster_first_pass) |> ungroup() |> arrange(desc(composite_score)) |>
    slice_head(n=quota) |> select(-dc)
  need <- quota-nrow(first)
  if(need>0) first <- bind_rows(first, eligible |> filter(!pathway %in% first$pathway) |> slice_head(n=need))
  first |> arrange(desc(composite_score)) |> mutate(selection_rank_within_block=row_number(), selected=TRUE)
}

selected <- stats |> group_by(functional_block) |> group_modify(~select_block(.x,.y$functional_block)) |> ungroup()
selected <- selected |> mutate(functional_block=factor(functional_block,levels=block_order)) |> arrange(functional_block,selection_rank_within_block)
stats_out <- stats |> left_join(selected |> select(pathway,selected,selection_rank_within_block),by="pathway") |>
  mutate(selected=replace_na(selected,FALSE),functional_block=factor(functional_block,levels=block_order)) |>
  arrange(functional_block,desc(selected),selection_rank_within_block,desc(composite_score))

write_csv(stats_out,file.path(output_dir,"04_all_pathway_selection_statistics_v3.csv"))
write_csv(selected,file.path(output_dir,"05_selected_representative_pathways_v3.csv"))
write_csv(selected |> count(functional_block,dominant_cluster,name="n_selected_pathways"),file.path(output_dir,"06_selected_pathways_dominant_cluster_summary_v3.csv"))

classification_for_supp <- classification |> select(pathway,pathway_label,collection_classification=collection,functional_block,eligible_for_v3_selection,exclusion_reason)
supp <- gsea |> left_join(classification_for_supp,by="pathway") |> mutate(
  collection=coalesce(as.character(collection),as.character(collection_classification),"Other"),
  pathway_label=coalesce(pathway_label,str_to_sentence(str_replace_all(pathway,"_"," "))),
  significant=padj_safe<fdr_cutoff & abs(NES)>=nes_cutoff
) |> select(functional_block,pathway,pathway_label,collection,eligible_for_v3_selection,exclusion_reason,cluster,NES,padj=padj_safe,significant)
write_csv(supp,file.path(output_dir,"07_Dataset_EV7_all_pathways_all_clusters_v3.csv"))
message("Selected approximately 40 neutrophil-focused representative pathways without hard-coding clusters.")
