clean_pathway_label <- function(x) {
  x |>
    stringr::str_remove("^HALLMARK_") |>
    stringr::str_remove("^REACTOME_") |>
    stringr::str_remove("^GOBP_") |>
    stringr::str_replace_all("_", " ") |>
    stringr::str_to_sentence()
}

row_order_by_nes <- function(mat) {
  if (nrow(mat) <= 2) return(rownames(mat))
  x <- mat
  x[!is.finite(x)] <- 0
  rownames(mat)[stats::hclust(stats::dist(x), method = "complete")$order]
}

prepare_version <- function(pathways, nes_complete, fdr_complete, cluster_order) {
  nes <- nes_complete[pathways, cluster_order, drop = FALSE]
  fdr <- fdr_complete[pathways, cluster_order, drop = FALSE]
  ord <- row_order_by_nes(nes)
  nes <- nes[ord, , drop = FALSE]
  fdr <- fdr[ord, , drop = FALSE]
  labels <- make.unique(vapply(rownames(nes), clean_pathway_label, character(1)), sep = " ")
  list(nes = nes, fdr = fdr, labels = labels)
}

significance_symbols <- function(nes, fdr, fdr_cutoff, nes_cutoff, require_nes = TRUE) {
  sig <- if (require_nes) fdr < fdr_cutoff & abs(nes) > nes_cutoff else fdr < fdr_cutoff
  out <- matrix("", nrow(nes), ncol(nes), dimnames = dimnames(nes))
  out[!is.na(sig) & sig] <- "\u25cf"
  out
}

save_heatmap <- function(version, name, title, output_dir, pathway_summary,
                         fdr_cutoff, nes_cutoff, nes_color_limit, require_nes) {
  nes <- version$nes
  fdr <- version$fdr
  symbols <- significance_symbols(nes, fdr, fdr_cutoff, nes_cutoff, require_nes)

  ann <- pathway_summary |>
    dplyr::filter(pathway %in% rownames(nes)) |>
    dplyr::select(pathway, collection) |>
    dplyr::distinct() |>
    tibble::column_to_rownames("pathway")
  ann <- ann[rownames(nes), , drop = FALSE]

  ann_colors <- list(collection = c(
    "Hallmark" = "#4C78A8", "Reactome" = "#F58518", "GO BP" = "#54A24B"
  ))
  breaks <- seq(-nes_color_limit, nes_color_limit, length.out = 101)
  colors <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
  n <- nrow(nes)
  pdf_h <- max(7, min(60, 2.8 + n * 0.18))
  png_h <- max(1600, 300 + n * 38)

  draw <- function() {
    pheatmap::pheatmap(
      nes, color = colors, breaks = breaks,
      cluster_rows = FALSE, cluster_cols = FALSE,
      display_numbers = symbols, number_color = "black", fontsize_number = 7,
      labels_row = version$labels, labels_col = paste0("Cluster ", colnames(nes)),
      annotation_row = ann, annotation_colors = ann_colors, annotation_names_row = FALSE,
      border_color = "grey85", fontsize = 9,
      fontsize_row = ifelse(n > 100, 5, ifelse(n > 70, 6, 7)),
      fontsize_col = 10, angle_col = 0, cellwidth = 58,
      cellheight = ifelse(n > 100, 8, ifelse(n > 70, 10, 13)),
      main = title, legend_breaks = c(-2, -1, 0, 1, 2), na_col = "grey90"
    )
  }

  grDevices::pdf(file.path(output_dir, paste0(name, ".pdf")), width = 12, height = pdf_h, useDingbats = FALSE)
  draw(); grDevices::dev.off()
  grDevices::png(file.path(output_dir, paste0(name, ".png")), width = 3000, height = png_h, res = 300)
  draw(); grDevices::dev.off()

  readr::write_csv(
    as.data.frame(nes) |> tibble::rownames_to_column("pathway") |>
      dplyr::mutate(pathway_label = version$labels, .after = pathway),
    file.path(output_dir, paste0(name, "_NES_matrix.csv"))
  )
  readr::write_csv(
    as.data.frame(fdr) |> tibble::rownames_to_column("pathway") |>
      dplyr::mutate(pathway_label = version$labels, .after = pathway),
    file.path(output_dir, paste0(name, "_FDR_matrix.csv"))
  )
}

save_network_plot_file <- function(publication_result, ranking, output_dir, label_top_n = 30) {
  g <- publication_result$graph
  if (igraph::vcount(g) == 0) return(invisible(NULL))

  assign <- publication_result$assignment
  rank_map <- setNames(ranking$priority_rank, ranking$pathway)
  igraph::V(g)$community <- assign$community[match(igraph::V(g)$name, assign$pathway)]
  igraph::V(g)$is_rep <- assign$retained[match(igraph::V(g)$name, assign$pathway)]
  igraph::V(g)$size <- ifelse(igraph::V(g)$is_rep, 7, 3)
  top_nodes <- names(sort(rank_map, decreasing = FALSE))[seq_len(min(label_top_n, length(rank_map)))]
  igraph::V(g)$label <- ifelse(
    igraph::V(g)$name %in% top_nodes,
    vapply(igraph::V(g)$name, clean_pathway_label, character(1)),
    NA
  )

  grDevices::pdf(file.path(output_dir, "GSEA_PUBLICATION_network.pdf"), width = 14, height = 12, useDingbats = FALSE)
  set.seed(1234)
  plot(
    g,
    vertex.label.cex = 0.55,
    vertex.frame.color = NA,
    vertex.color = as.integer(as.factor(igraph::V(g)$community)),
    edge.width = 0.5 + 2 * igraph::E(g)$weight,
    edge.color = grDevices::adjustcolor("grey50", alpha.f = 0.35),
    layout = igraph::layout_with_fr(g, weights = igraph::E(g)$weight)
  )
  grDevices::dev.off()
}
