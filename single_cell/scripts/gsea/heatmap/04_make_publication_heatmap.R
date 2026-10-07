# Manuscript output: Fig. EV4E (broader functional GSEA landscape)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(scales)
})
source("config.R")

selected <- read_csv(
  file.path(output_dir, "05_selected_representative_pathways_v3.csv"),
  show_col_types = FALSE
)

gsea <- read_csv(
  file.path(output_dir, "01_all_clusters_GSEA_long.csv"),
  show_col_types = FALSE
) |>
  mutate(
    cluster = factor(as.character(cluster), levels = cluster_order),
    NES = as.numeric(NES),
    padj = pmax(as.numeric(padj), 1e-300)
  )

plot_df <- selected |>
  select(pathway, pathway_label, functional_block, selection_rank_within_block) |>
  left_join(gsea |> select(pathway, cluster, NES, padj), by = "pathway") |>
  mutate(
    pathway_label = if_else(
      is.na(pathway_label) | pathway_label == "",
      str_to_sentence(str_replace_all(pathway, "_", " ")),
      pathway_label
    ),
    significant = padj < fdr_cutoff &
      if (mark_requires_nes_cutoff) abs(NES) >= nes_cutoff else TRUE
  )

# Return pathway order for one functional block.
# Clustering changes only the row display order; NES and FDR values are unchanged.
order_pathways_in_block <- function(df,
                                    cluster_levels,
                                    use_clustering = TRUE,
                                    distance_method = "spearman",
                                    linkage_method = "complete") {
  meta <- df |>
    distinct(pathway, pathway_label, selection_rank_within_block) |>
    arrange(selection_rank_within_block, pathway_label)

  if (!use_clustering || nrow(meta) <= 2L) {
    return(meta |> mutate(clustered_order_within_block = row_number()))
  }

  mat_df <- df |>
    mutate(cluster = factor(as.character(cluster), levels = cluster_levels)) |>
    select(pathway, cluster, NES) |>
    group_by(pathway, cluster) |>
    summarise(
      NES = {
        z <- NES[is.finite(NES)]
        if (length(z) == 0L) NA_real_ else z[which.max(abs(z))]
      },
      .groups = "drop"
    ) |>
    pivot_wider(
      names_from = cluster,
      values_from = NES,
      values_fill = NA_real_
    ) |>
    right_join(meta |> select(pathway), by = "pathway")

  missing_cluster_cols <- setdiff(cluster_levels, colnames(mat_df))
  if (length(missing_cluster_cols) > 0L) {
    for (nm in missing_cluster_cols) mat_df[[nm]] <- NA_real_
  }

  mat <- as.matrix(mat_df[, cluster_levels, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- mat_df$pathway
  mat[!is.finite(mat)] <- 0

  # Remove rows with no variation from distance calculation only.
  # They are appended in their original selection order.
  variable <- apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
  variable[is.na(variable)] <- FALSE

  ordered_paths <- character(0)
  if (sum(variable) >= 2L) {
    mat_var <- mat[variable, , drop = FALSE]

    if (identical(distance_method, "spearman")) {
      cor_mat <- suppressWarnings(stats::cor(t(mat_var), method = "spearman", use = "pairwise.complete.obs"))
      cor_mat[!is.finite(cor_mat)] <- 0
      diag(cor_mat) <- 1
      # Preserve the square matrix dimensions. Using pmax() directly can
      # drop the dim attribute and produce an invalid dist object.
      dist_mat <- 1 - cor_mat
      dist_mat[!is.finite(dist_mat)] <- 1
      dist_mat[dist_mat < 0] <- 0
      diag(dist_mat) <- 0
      d <- stats::as.dist(dist_mat)
    } else if (identical(distance_method, "euclidean")) {
      d <- stats::dist(mat_var, method = "euclidean")
    } else {
      stop("within_block_distance must be either 'spearman' or 'euclidean'.")
    }

    expected_d_length <- nrow(mat_var) * (nrow(mat_var) - 1L) / 2L
    if (length(d) != expected_d_length || any(!is.finite(d))) {
      warning(
        "Invalid within-block distance object; using original selection order for this block."
      )
      ordered_paths <- rownames(mat_var)
    } else {
      hc <- stats::hclust(d, method = linkage_method)
      ordered_paths <- rownames(mat_var)[hc$order]
    }
  } else if (sum(variable) == 1L) {
    ordered_paths <- rownames(mat)[variable]
  }

  nonvariable_paths <- meta |>
    filter(!pathway %in% ordered_paths) |>
    arrange(selection_rank_within_block, pathway_label) |>
    pull(pathway)

  final_paths <- c(ordered_paths, nonvariable_paths)

  tibble(
    pathway = final_paths,
    clustered_order_within_block = seq_along(final_paths)
  ) |>
    left_join(meta, by = "pathway")
}

# Normalize labels only for robust matching of the fixed manual order.
normalize_label <- function(x) {
  x |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", " ") |>
    stringr::str_squish()
}

order_pathways_manual <- function(df, requested_labels) {
  meta <- df |>
    distinct(pathway, pathway_label, selection_rank_within_block) |>
    arrange(selection_rank_within_block, pathway_label) |>
    mutate(label_key = normalize_label(pathway_label))

  requested_keys <- normalize_label(requested_labels)
  matched <- character(0)
  for (key in requested_keys) {
    hit <- meta |> filter(label_key == key) |> pull(pathway)
    if (length(hit) > 0L) matched <- c(matched, hit[[1]])
  }

  # Any pathway not listed explicitly is appended in its original selection order.
  remaining <- meta |>
    filter(!pathway %in% matched) |>
    arrange(selection_rank_within_block, pathway_label) |>
    pull(pathway)

  final_paths <- c(matched, remaining)
  tibble(
    pathway = final_paths,
    clustered_order_within_block = seq_along(final_paths),
    ordering_method = if_else(final_paths %in% matched, "manual_fixed", "manual_block_unlisted_append")
  ) |>
    left_join(meta |> select(-label_key), by = "pathway")
}

block_order_table <- plot_df |>
  mutate(functional_block = factor(functional_block, levels = block_order)) |>
  group_by(functional_block) |>
  group_modify(~ {
    block_name <- as.character(.y$functional_block[[1]])
    if (block_name %in% manual_order_blocks) {
      order_pathways_manual(.x, manual_pathway_order[[block_name]])
    } else {
      order_pathways_in_block(
        .x,
        cluster_levels = cluster_order,
        use_clustering = cluster_pathways_within_block,
        distance_method = within_block_distance,
        linkage_method = within_block_linkage
      ) |>
        mutate(ordering_method = if_else(
          cluster_pathways_within_block,
          paste0("hierarchical_", within_block_distance),
          "selection_order"
        ))
    }
  }) |>
  ungroup() |>
  mutate(functional_block = factor(functional_block, levels = block_order)) |>
  arrange(functional_block, clustered_order_within_block)

# ggplot draws the first factor level at the bottom, so reverse the complete
# top-to-bottom pathway sequence while preserving the block order.
pathway_levels_top_to_bottom <- block_order_table |> pull(pathway_label)
pathway_levels_for_ggplot <- rev(pathway_levels_top_to_bottom)

plot_df <- plot_df |>
  left_join(
    block_order_table |>
      select(pathway, functional_block, clustered_order_within_block),
    by = c("pathway", "functional_block")
  ) |>
  mutate(
    functional_block = factor(functional_block, levels = block_order),
    pathway_label = factor(pathway_label, levels = pathway_levels_for_ggplot)
  )

p <- ggplot(
  plot_df,
  aes(
    cluster,
    pathway_label,
    fill = pmax(pmin(NES, nes_color_limit), -nes_color_limit)
  )
) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_point(
    data = filter(plot_df, significant),
    shape = 21,
    size = 1.65,
    stroke = 0.35,
    fill = "white"
  ) +
  facet_grid(
    functional_block ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y",
    drop = TRUE
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-nes_color_limit, nes_color_limit),
    oob = scales::squish,
    name = "NES"
  ) +
  labs(
    x = "Cluster",
    y = NULL,
    title = "Neutrophil-focused GSEA programs across clusters",
    subtitle = paste0(
      "White dot: FDR < ", fdr_cutoff,
      " and |NES| >= ", nes_cutoff,
      "; Cell cycle and IFN ordered manually; other blocks clustered within block (",
      within_block_distance, ")"
    )
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey96", color = "grey75"),
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

height <- max(
  heatmap_min_height,
  n_distinct(plot_df$pathway) * heatmap_height_per_pathway
)

ggsave(
  file.path(output_dir, "11_publication_heatmap_all_clusters_v3_3.pdf"),
  p,
  width = heatmap_width,
  height = height,
  device = cairo_pdf
)
ggsave(
  file.path(output_dir, "11_publication_heatmap_all_clusters_v3_3.png"),
  p,
  width = heatmap_width,
  height = height,
  dpi = 400
)

write_csv(
  plot_df,
  file.path(output_dir, "11_publication_heatmap_plot_data_v3_3.csv")
)
write_csv(
  block_order_table |>
    mutate(
      clustering_enabled = cluster_pathways_within_block,
      distance_method = within_block_distance,
      linkage_method = within_block_linkage
    ),
  file.path(output_dir, "11_pathway_order_within_blocks_v3_3.csv")
)

message("v3.3 publication heatmap saved: Cell cycle and Interferon fixed manually; other blocks clustered within block.")

