read_gsea_file <- function(file_path, cluster_id) {
  x <- readr::read_csv(file_path, show_col_types = FALSE)

  required <- c("pathway", "pval", "padj", "NES", "size", "leadingEdge", "collection")
  missing <- setdiff(required, colnames(x))
  if (length(missing) > 0) {
    stop("Missing columns in ", basename(file_path), ": ", paste(missing, collapse = ", "))
  }

  x |>
    dplyr::transmute(
      cluster = as.character(cluster_id),
      pathway = as.character(pathway),
      collection_raw = as.character(collection),
      pval = as.numeric(pval),
      padj = as.numeric(padj),
      NES = as.numeric(NES),
      size = as.numeric(size),
      leadingEdge = as.character(leadingEdge)
    )
}

read_all_gsea <- function(input_files, cluster_order) {
  missing_files <- input_files[!file.exists(input_files)]
  if (length(missing_files) > 0) {
    stop("Input files not found:\n", paste(missing_files, collapse = "\n"))
  }

  x <- purrr::imap_dfr(input_files, ~ read_gsea_file(.x, .y)) |>
    dplyr::mutate(
      collection = dplyr::case_when(
        stringr::str_detect(pathway, "^HALLMARK_") ~ "Hallmark",
        stringr::str_detect(pathway, "^REACTOME_") ~ "Reactome",
        stringr::str_detect(pathway, "^GOBP_") ~ "GO BP",
        stringr::str_detect(stringr::str_to_lower(collection_raw), "hallmark") ~ "Hallmark",
        stringr::str_detect(stringr::str_to_lower(collection_raw), "reactome") ~ "Reactome",
        stringr::str_detect(stringr::str_to_lower(collection_raw), "go") ~ "GO BP",
        TRUE ~ collection_raw
      ),
      cluster = factor(cluster, levels = cluster_order)
    )

  dup <- x |>
    dplyr::count(cluster, pathway) |>
    dplyr::filter(n > 1)

  if (nrow(dup) > 0) {
    warning("Duplicate pathway-cluster rows found; retaining the lowest-FDR row.")
    x <- x |>
      dplyr::arrange(cluster, pathway, padj, dplyr::desc(abs(NES))) |>
      dplyr::distinct(cluster, pathway, .keep_all = TRUE)
  }

  x
}

summarize_pathways <- function(gsea_all, fdr_cutoff, nes_cutoff) {
  gsea_all |>
    dplyr::group_by(pathway) |>
    dplyr::summarise(
      collection = dplyr::first(collection),
      min_FDR = ifelse(all(is.na(padj)), NA_real_, min(padj, na.rm = TRUE)),
      max_abs_NES = ifelse(all(is.na(NES)), NA_real_, max(abs(NES), na.rm = TRUE)),
      representative_size = {
        idx <- which(is.finite(padj))
        if (length(idx) == 0) dplyr::first(size) else size[idx[which.min(padj[idx])]]
      },
      n_clusters_full_cutoff = sum(padj < fdr_cutoff & abs(NES) > nes_cutoff, na.rm = TRUE),
      include_complete = any(padj < fdr_cutoff & abs(NES) > nes_cutoff, na.rm = TRUE),
      .groups = "drop"
    )
}
