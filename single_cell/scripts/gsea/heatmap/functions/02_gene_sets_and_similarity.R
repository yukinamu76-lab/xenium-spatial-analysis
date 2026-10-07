safe_msigdbr <- function(collection, subcollection = NULL) {
  tryCatch({
    if (is.null(subcollection)) {
      msigdbr::msigdbr(db_species = "MM", species = "Mus musculus", collection = collection)
    } else {
      msigdbr::msigdbr(
        db_species = "MM", species = "Mus musculus",
        collection = collection, subcollection = subcollection
      )
    }
  }, error = function(e) {
    warning("msigdbr retrieval failed for ", collection, ": ", conditionMessage(e))
    tibble::tibble()
  })
}

parse_leading_edge <- function(x) {
  if (is.na(x) || x == "") return(character())
  x <- stringr::str_remove_all(x, "^c\\(|\\)$")
  genes <- stringr::str_split(x, "[,;/|[:space:]]+")[[1]]
  genes <- stringr::str_remove_all(genes, "[\"']")
  unique(genes[!is.na(genes) & genes != ""])
}

build_pathway_gene_sets <- function(gsea_all, pathways) {
  
  hallmark <- safe_msigdbr("MH")
  reactome <- safe_msigdbr("M2", "CP:REACTOME")
  gobp <- safe_msigdbr("M5", "GO:BP")
  
  msig_all <- dplyr::bind_rows(
    hallmark,
    reactome,
    gobp
  )
  
  gene_col <- if ("gene_symbol" %in% names(msig_all)) {
    
    "gene_symbol"
    
  } else if ("human_gene_symbol" %in% names(msig_all)) {
    
    "human_gene_symbol"
    
  } else {
    
    NA_character_
  }
  
  msig_sets <- list()
  
  if (
    nrow(msig_all) > 0 &&
    !is.na(gene_col) &&
    "gs_name" %in% names(msig_all)
  ) {
    
    msig_sets_df <- msig_all |>
      dplyr::filter(
        gs_name %in% pathways
      ) |>
      dplyr::select(
        gs_name,
        dplyr::all_of(gene_col)
      ) |>
      dplyr::filter(
        !is.na(.data[[gene_col]]),
        .data[[gene_col]] != ""
      ) |>
      dplyr::distinct()
    
    msig_sets <- split(
      msig_sets_df[[gene_col]],
      msig_sets_df$gs_name
    ) |>
      purrr::map(unique)
  }
  
  le_sets <- gsea_all |>
    dplyr::filter(
      pathway %in% pathways
    ) |>
    dplyr::group_by(pathway) |>
    dplyr::summarise(
      genes = list(
        unique(
          unlist(
            purrr::map(
              leadingEdge,
              parse_leading_edge
            )
          )
        )
      ),
      .groups = "drop"
    ) |>
    tibble::deframe()
  
  sets <- stats::setNames(
    purrr::map(
      pathways,
      function(pw) {
        
        x <- msig_sets[[pw]]
        
        if (
          !is.null(x) &&
          length(x) > 0
        ) {
          
          unique(x)
          
        } else {
          
          unique(le_sets[[pw]])
        }
      }
    ),
    pathways
  )
  
  source <- tibble::tibble(
    pathway = pathways,
    
    gene_set_source = purrr::map_chr(
      pathways,
      function(pw) {
        
        x <- msig_sets[[pw]]
        
        if (
          !is.null(x) &&
          length(x) > 0
        ) {
          
          "Full MSigDB gene set"
          
        } else {
          
          "Leading-edge fallback"
        }
      }
    ),
    
    n_genes_used = purrr::map_int(
      sets,
      length
    )
  )
  
  list(
    sets = sets,
    source = source
  )
}

jaccard_similarity <- function(a, b) {
  a <- unique(a[a != "" & !is.na(a)])
  b <- unique(b[b != "" & !is.na(b)])
  u <- length(union(a, b))
  if (u == 0) return(NA_real_)
  length(intersect(a, b)) / u
}

overlap_coefficient <- function(a, b) {
  a <- unique(a[a != "" & !is.na(a)])
  b <- unique(b[b != "" & !is.na(b)])
  d <- min(length(a), length(b))
  if (d == 0) return(NA_real_)
  length(intersect(a, b)) / d
}

make_numeric_matrix <- function(data, value_column, pathways, cluster_order) {
  wide <- data |>
    dplyr::filter(pathway %in% pathways) |>
    dplyr::mutate(cluster = as.character(cluster)) |>
    dplyr::select(pathway, cluster, dplyr::all_of(value_column)) |>
    tidyr::pivot_wider(names_from = cluster, values_from = dplyr::all_of(value_column))

  for (cl in setdiff(cluster_order, names(wide))) wide[[cl]] <- NA_real_

  mat <- wide |>
    dplyr::select(pathway, dplyr::all_of(cluster_order)) |>
    tibble::column_to_rownames("pathway") |>
    as.matrix()
  storage.mode(mat) <- "numeric"
  mat[pathways, cluster_order, drop = FALSE]
}

calculate_pairwise_similarity <- function(pathways, gene_sets, nes_matrix, cluster_order) {
  if (length(pathways) < 2) {
    return(tibble::tibble(
      pathway_A = character(), pathway_B = character(),
      jaccard = numeric(), overlap_coefficient = numeric(), NES_spearman = numeric()
    ))
  }

  pairs <- combn(pathways, 2, simplify = FALSE)
  purrr::map_dfr(pairs, function(pair) {
    a <- pair[1]; b <- pair[2]
    va <- as.numeric(nes_matrix[a, cluster_order])
    vb <- as.numeric(nes_matrix[b, cluster_order])
    ok <- is.finite(va) & is.finite(vb)
    rho <- if (sum(ok) >= 3 && length(unique(va[ok])) > 1 && length(unique(vb[ok])) > 1) {
      suppressWarnings(stats::cor(va[ok], vb[ok], method = "spearman"))
    } else NA_real_

    tibble::tibble(
      pathway_A = a,
      pathway_B = b,
      jaccard = jaccard_similarity(gene_sets[[a]], gene_sets[[b]]),
      overlap_coefficient = overlap_coefficient(gene_sets[[a]], gene_sets[[b]]),
      NES_spearman = rho
    )
  })
}
