## ============================================================
## Fig. 6B pooled cluster-composition summary
##
## This public-release reproduction script implements the documented
## calculation used for Fig. 6B: pool cells across the two biological
## replicates within each treatment group, calculate the fraction of
## neutrophils in each cluster, and report Alb/PBS for each cluster.
## ============================================================
source(file.path("scripts", "scRNA", "00_setup.R"))

neut <- readRDS(file.path(dir_rds, "07_neutrophil_clustered.rds"))
out_dir <- file.path(dir_output, "07_neutrophil_composition")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta <- neut@meta.data |>
  tibble::as_tibble(rownames = "cell") |>
  dplyr::mutate(cluster = as.character(seurat_clusters), group = as.character(group))

pooled <- meta |>
  dplyr::count(group, cluster, name = "n_cells") |>
  dplyr::group_by(group) |>
  dplyr::mutate(total_neutrophils = sum(n_cells), proportion = n_cells / total_neutrophils) |>
  dplyr::ungroup()

wide <- pooled |>
  dplyr::select(group, cluster, proportion) |>
  tidyr::pivot_wider(names_from = group, values_from = proportion, values_fill = 0)

if (!all(c("PBS", "Alb") %in% colnames(wide))) {
  stop("Expected group labels PBS and Alb were not found.")
}

fig6b <- wide |>
  dplyr::mutate(Alb_over_PBS = Alb / PBS) |>
  dplyr::arrange(as.integer(cluster))

readr::write_csv(pooled, file.path(out_dir, "Fig6B_pooled_cluster_proportions.csv"))
readr::write_csv(fig6b, file.path(out_dir, "Fig6B_Alb_over_PBS_cluster_ratio.csv"))
