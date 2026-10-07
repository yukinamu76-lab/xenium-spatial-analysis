suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(stringr); library(tibble)
})
source("config.R")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path("functions", "01_read_data.R"))

gsea_all <- read_all_gsea(input_files, cluster_order) |>
  mutate(cluster = factor(as.character(cluster), levels = cluster_order))
write_csv(gsea_all, file.path(output_dir, "01_all_clusters_GSEA_long.csv"))

matches_any <- function(x, patterns) str_detect(str_to_upper(x), paste(patterns, collapse = "|"))
detect_collection <- function(pathway) case_when(
  str_starts(pathway, "HALLMARK_") ~ "Hallmark",
  str_starts(pathway, "REACTOME_") ~ "Reactome",
  str_starts(pathway, "GOBP_") ~ "GO BP",
  TRUE ~ "Other"
)

assign_block_v3 <- function(pathway) {
  x <- str_to_upper(pathway)
  if (matches_any(x, c("CELL_CYCLE", "MITOTIC", "M_PHASE", "G2M", "E2F", "CHROMOSOME", "CHROMATID", "CELL_DIVISION", "DNA_REPLICATION", "SPINDLE", "CYTOKINESIS", "KINETOCHORE", "CENTROSOME"))) return("Cell cycle")
  if (matches_any(x, c("INTERFERON", "TYPE_I_INTERFERON", "TYPE_II_INTERFERON", "ANTIVIRAL", "RESPONSE_TO_VIRUS", "VIRAL_PROCESS"))) return("Interferon response")
  if (matches_any(x, c("CHEMOTAX", "LEUKOCYTE_MIGRATION", "CELL_MIGRATION", "MOTILITY", "ADHESION", "INTEGRIN", "EXTRACELLULAR_MATRIX", "ECM"))) return("Migration and chemotaxis")
  if (matches_any(x, c("INFLAM", "TOLL_LIKE", "NOD_LIKE", "NLR", "DEFENSE_RESPONSE", "IMMUNE_RESPONSE", "BACTERIAL", "ANTIMICROBIAL", "PHAGOCYT", "CYTOKINE", "LEUKOCYTE_ACTIVATION", "DEGRANULATION", "GRANULE", "EXOCYTOSIS", "RESPIRATORY_BURST", "COMPLEMENT"))) return("Innate immunity")
  if (matches_any(x, c("TRANSLATION", "RIBOSOM", "RRNA", "MRNA_PROCESS", "RNA_PROCESS", "RNA_SPLIC", "SPLICEOSOME", "RNA_STABILITY", "RNA_CATABOLIC"))) return("RNA processing and translation")
  if (matches_any(x, c("PROTEIN_FOLDING", "UNFOLDED_PROTEIN", "PROTEOSTASIS", "PROTEASOME", "PROTEIN_DEGRADATION", "PROTEIN_CATABOLIC", "AUTOPHAG", "RESPONSE_TO_STRESS", "OXIDATIVE_STRESS", "ER_STRESS", "HEAT_SHOCK"))) return("Cellular stress and protein homeostasis")
  if (matches_any(x, c("OXIDATIVE_PHOSPHORYLATION", "ELECTRON_TRANSPORT", "MITOCHONDR", "GLYCOLYSIS", "FATTY_ACID", "CHOLESTEROL", "LIPID_METAB", "AMINO_ACID", "CARBOHYDRATE", "TCA_CYCLE", "PENTOSE", "HEME_METAB", "METABOL"))) return("Metabolism")
  if (matches_any(x, c("SIGNAL_TRANSDUCTION", "MAPK", "PI3K", "AKT", "MTOR", "JAK_STAT", "NF_KAPPA", "P53", "PTEN", "KRAS", "FOXO", "APOPT", "CELL_DEATH", "SURVIVAL", "DIFFERENTIATION", "PROLIFERATION"))) return("Signal transduction and cell fate")
  NA_character_
}

# Manual classifications are retained only when they map cleanly to the v3 blocks.
manual <- read_csv(classification_file, show_col_types = FALSE) |>
  distinct(pathway, .keep_all = TRUE) |>
  mutate(
    pathway_label = coalesce(as.character(pathway_label), str_to_sentence(str_replace_all(pathway, "_", " "))),
    collection = coalesce(as.character(collection), detect_collection(pathway)),
    v3_block_from_name = vapply(pathway, assign_block_v3, character(1)),
    functional_block_v3 = v3_block_from_name,
    classification_source = "v3 biological keyword classification"
  ) |>
  select(-v3_block_from_name)

all_pathways <- gsea_all |> distinct(pathway)
classification <- all_pathways |>
  left_join(manual |> select(pathway, pathway_label, collection, manual_note), by = "pathway") |>
  mutate(
    pathway_label = coalesce(pathway_label, str_to_sentence(str_replace_all(pathway, "_", " "))),
    collection = coalesce(collection, vapply(pathway, detect_collection, character(1))),
    functional_block = vapply(pathway, assign_block_v3, character(1)),
    excluded_tissue_specific = matches_any(pathway, exclude_patterns),
    excluded_generic = matches_any(pathway, generic_patterns),
    classification_source = "v3 biological keyword classification",
    manual_note = coalesce(manual_note, "")
  ) |>
  mutate(
    exclusion_reason = case_when(
      excluded_tissue_specific ~ "Tissue-specific/developmental term",
      excluded_generic ~ "Too generic for compact publication heatmap",
      is.na(functional_block) ~ "Not assigned to a neutrophil-focused v3 block",
      TRUE ~ ""
    ),
    eligible_for_v3_selection = exclusion_reason == ""
  )

write_csv(classification, file.path(output_dir, "02_classification_used_v3.csv"))
write_csv(classification |> filter(!eligible_for_v3_selection), file.path(output_dir, "02b_excluded_or_unassigned_pathways_for_review.csv"))
write_csv(classification |> filter(eligible_for_v3_selection), file.path(output_dir, "02c_v3_eligible_pathways.csv"))

wide_nes <- gsea_all |> select(pathway, cluster, NES) |> distinct() |> pivot_wider(names_from = cluster, values_from = NES, names_prefix = "NES_C")
wide_fdr <- gsea_all |> select(pathway, cluster, padj) |> distinct() |> pivot_wider(names_from = cluster, values_from = padj, names_prefix = "FDR_C")
full_table <- classification |> left_join(wide_nes, by = "pathway") |> left_join(wide_fdr, by = "pathway")
write_csv(full_table, file.path(output_dir, "03_all_pathways_NES_FDR_wide_v3.csv"))
message("Built v3 all-cluster tables and neutrophil-focused classification.")
