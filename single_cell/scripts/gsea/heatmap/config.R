# ============================================================
# GSEA all-cluster pipeline v3: neutrophil-focused publication version
# ============================================================
input_dir <- file.path(getwd(), "Outputs", "10_GSEA")
output_dir <- file.path(input_dir, "GSEA_allClusters_v3_output")
classification_file <- file.path("classification", "STEP1_pathway_classification.csv")
cluster_order <- as.character(0:8)
input_files <- stats::setNames(
  file.path(input_dir, paste0("cluster", cluster_order), paste0("cluster", cluster_order, "_GSEA_all_collections_results.csv")),
  cluster_order
)

fdr_cutoff <- 0.05
nes_cutoff <- 1.5

# Publication-focused blocks and quotas: total target = 40 pathways.
block_order <- c(
  "Cell cycle",
  "Interferon response",
  "Innate immunity",
  "Migration and chemotaxis",
  "Metabolism",
  "RNA processing and translation",
  "Cellular stress and protein homeostasis",
  "Signal transduction and cell fate"
)
pathways_per_block <- c(
  "Cell cycle" = 6,
  "Interferon response" = 4,
  "Innate immunity" = 6,
  "Migration and chemotaxis" = 5,
  "Metabolism" = 5,
  "RNA processing and translation" = 5,
  "Cellular stress and protein homeostasis" = 4,
  "Signal transduction and cell fate" = 5
)
default_pathways_per_block <- 5

selection_weights <- c(
  significance = 0.30,
  effect = 0.30,
  specificity = 0.25,
  rarity = 0.10,
  interpretability = 0.05
)
max_per_dominant_cluster_first_pass <- 1

# Exclude tissue-specific or poorly interpretable terms from representative selection.
# They remain in the complete results reported in Dataset EV7.
exclude_patterns <- c(
  "CARDIAC", "HEART", "MUSCLE_CELL", "SKELETAL_MUSCLE", "SMOOTH_MUSCLE",
  "MEGAKARYOCYTE", "PLATELET", "ERYTHROCYTE", "ERYTHROPOIESIS",
  "EPITHELIUM", "EPITHELIAL", "KERATIN", "SKIN", "HAIR", "BONE", "CARTILAGE",
  "NEURON", "NEURAL", "AXON", "SYNAP", "RETINA", "EYE", "EAR",
  "KIDNEY", "RENAL", "LIVER_DEVELOPMENT", "PANCREAS", "LUNG_DEVELOPMENT",
  "REPRODUCT", "GONAD", "PLACENTA", "EMBRYONIC", "MORPHOGENESIS",
  "ANATOMICAL_STRUCTURE", "ORGAN_DEVELOPMENT", "TISSUE_DEVELOPMENT"
)

# Terms that are usually too generic to be informative in a compact heatmap.
generic_patterns <- c(
  "REGULATION_OF_BIOLOGICAL_PROCESS$", "CELLULAR_PROCESS$", "METABOLIC_PROCESS$",
  "DEVELOPMENTAL_PROCESS$", "RESPONSE_TO_STIMULUS$", "SIGNALING$",
  "POSITIVE_REGULATION_OF_CELLULAR_PROCESS$", "NEGATIVE_REGULATION_OF_CELLULAR_PROCESS$"
)


# Pathway ordering within each functional block.
# TRUE: hierarchically cluster pathways using their NES patterns across clusters.
cluster_pathways_within_block <- TRUE
# "spearman" groups pathways with similar NES shapes, regardless of absolute magnitude.
# "euclidean" groups pathways using absolute NES differences.
within_block_distance <- "spearman"
within_block_linkage <- "complete"

# Hybrid pathway ordering for the publication heatmap.
# Cell cycle and Interferon response use a fixed biological order;
# all other blocks retain within-block hierarchical clustering.
manual_order_blocks <- c("Cell cycle", "Interferon response")
manual_pathway_order <- list(
  "Cell cycle" = c(
    "E2f targets",
    "Gobp chromosome organization involved in meiotic cell cycle",
    "Gobp chromosome localization",
    "Reactome recruitment of mitotic centrosome proteins and complexes",
    "Gobp positive regulation of chromosome separation",
    "Mitotic cytokinesis"
  ),
  "Interferon response" = c(
    "Hallmark interferon alpha response",
    "Gobp positive regulation of type i interferon production",
    "Reactome ddx58 ifih1 mediated induction of interferon alpha beta",
    "Interferon gamma response"
  )
)

# Heatmap settings.
nes_color_limit <- 2.5
mark_requires_nes_cutoff <- TRUE
heatmap_width <- 12.5
heatmap_height_per_pathway <- 0.28
heatmap_min_height <- 8.5

reference_program_patterns <- list(
  Cell_cycle = "CELL_CYCLE|MITOTIC|CHROMOSOME|DNA_REPLICATION|SPINDLE|CYTOKINESIS|E2F|G2M",
  Interferon = "INTERFERON|IFN_ALPHA|IFN_GAMMA|ANTIVIRAL|RESPONSE_TO_VIRUS",
  Chemotaxis = "CHEMOTAX|MIGRATION|LEUKOCYTE_MIGRATION|CELL_ADHESION|INTEGRIN",
  Respiratory_burst = "RESPIRATORY_BURST|REACTIVE_OXYGEN|SUPEROXIDE|NADPH_OXIDASE",
  Translation = "TRANSLATION|RIBOSOM|RRNA_PROCESS|PROTEIN_SYNTHESIS",
  RNA_processing = "RNA_SPLIC|MRNA_PROCESS|SPLICEOSOME|RNA_PROCESSING",
  OXPHOS = "OXIDATIVE_PHOSPHORYLATION|RESPIRATORY_ELECTRON_TRANSPORT|MITOCHONDRIAL_RESPIRATION",
  Degranulation = "DEGRANULATION|GRANULE|EXOCYTOSIS",
  Apoptosis = "APOPTOTIC|APOPTOSIS|CELL_DEATH"
)
