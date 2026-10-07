make_representative_ranking <- function(complete_summary) {
  collection_rank <- c("Hallmark" = 1, "Reactome" = 2, "GO BP" = 3)
  complete_summary |>
    dplyr::mutate(
      collection_rank = unname(collection_rank[collection]),
      collection_rank = tidyr::replace_na(collection_rank, 99)
    ) |>
    dplyr::arrange(min_FDR, dplyr::desc(max_abs_NES), representative_size, collection_rank, pathway) |>
    dplyr::mutate(priority_rank = dplyr::row_number())
}

reduce_intermediate <- function(pairwise, ranking, jaccard_cutoff, overlap_cutoff, nes_cor_cutoff) {
  lookup <- pairwise |>
    dplyr::filter(
      ((is.finite(jaccard) & jaccard >= jaccard_cutoff) |
       (is.finite(overlap_coefficient) & overlap_coefficient >= overlap_cutoff)) &
      is.finite(NES_spearman) & NES_spearman >= nes_cor_cutoff
    )

  is_pair <- function(a, b) {
    any((lookup$pathway_A == a & lookup$pathway_B == b) |
        (lookup$pathway_A == b & lookup$pathway_B == a))
  }

  retained <- character()
  rows <- list()

  for (pw in ranking$pathway) {
    matches <- retained[purrr::map_lgl(retained, ~ is_pair(pw, .x))]
    if (length(matches) == 0) {
      retained <- c(retained, pw)
      rep_pw <- pw
      keep <- TRUE
    } else {
      rep_pw <- matches[1]
      keep <- FALSE
    }
    rows[[length(rows) + 1]] <- tibble::tibble(
      pathway = pw,
      representative_pathway = rep_pw,
      retained = keep
    )
  }

  assignment <- dplyr::bind_rows(rows) |>
    dplyr::left_join(ranking, by = "pathway")

  list(
    pathways = assignment |> dplyr::filter(retained) |> dplyr::pull(pathway),
    assignment = assignment,
    pairs = lookup
  )
}

build_publication_edge_table <- function(pairwise, candidate_pathways,
                                         min_nes_cor, jaccard_floor, overlap_floor,
                                         gene_weight, nes_weight) {
  pairwise |>
    dplyr::filter(pathway_A %in% candidate_pathways, pathway_B %in% candidate_pathways) |>
    dplyr::mutate(
      gene_similarity = pmax(
        dplyr::if_else(is.finite(jaccard), jaccard, 0),
        dplyr::if_else(is.finite(overlap_coefficient), overlap_coefficient, 0)
      ),
      nes_similarity = dplyr::if_else(
        is.finite(NES_spearman), pmax(0, NES_spearman), 0
      ),
      edge_score = gene_weight * gene_similarity + nes_weight * nes_similarity
    ) |>
    dplyr::filter(
      is.finite(NES_spearman), NES_spearman >= min_nes_cor,
      (jaccard >= jaccard_floor | overlap_coefficient >= overlap_floor)
    )
}

community_from_threshold <- function(nodes, edge_table, threshold, ranking, seed = 1234) {
  edges <- edge_table |>
    dplyr::filter(edge_score >= threshold) |>
    dplyr::select(pathway_A, pathway_B, edge_score)

  g <- igraph::make_empty_graph(n = 0, directed = FALSE) |>
    igraph::add_vertices(length(nodes), name = nodes)

  if (nrow(edges) > 0) {
    edge_vec <- as.vector(t(as.matrix(edges[, c("pathway_A", "pathway_B")])))
    g <- igraph::add_edges(g, edge_vec)
    igraph::E(g)$weight <- edges$edge_score
  }

  membership <- stats::setNames(as.character(seq_along(nodes)), nodes)

  if (igraph::ecount(g) > 0) {
    comp <- igraph::components(g)$membership
    for (cc in unique(comp)) {
      verts <- names(comp)[comp == cc]
      if (length(verts) == 1) next
      sg <- igraph::induced_subgraph(g, verts)
      set.seed(seed)
      cl <- igraph::cluster_louvain(sg, weights = igraph::E(sg)$weight)
      membership[verts] <- paste0("C", cc, "_", igraph::membership(cl))
    }
  }

  groups <- tibble::tibble(pathway = nodes, community = unname(membership[nodes])) |>
    dplyr::left_join(ranking, by = "pathway")

  reps <- groups |>
    dplyr::arrange(priority_rank) |>
    dplyr::group_by(community) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::transmute(community, representative_pathway = pathway)

  assignment <- groups |>
    dplyr::left_join(reps, by = "community") |>
    dplyr::mutate(retained = pathway == representative_pathway)

  list(
    threshold = threshold,
    graph = g,
    edges = edges,
    assignment = assignment,
    n_representatives = sum(assignment$retained)
  )
}

choose_publication_network <- function(nodes, edge_table, thresholds, ranking,
                                       target_count, count_range, seed = 1234) {
  trials <- purrr::map(thresholds, ~ community_from_threshold(
    nodes, edge_table, .x, ranking, seed
  ))

  summary <- purrr::map_dfr(trials, ~ tibble::tibble(
    threshold = .x$threshold,
    n_representatives = .x$n_representatives,
    n_edges = nrow(.x$edges)
  )) |>
    dplyr::mutate(
      in_target_range = n_representatives >= count_range[1] & n_representatives <= count_range[2],
      distance_to_target = abs(n_representatives - target_count)
    ) |>
    dplyr::arrange(dplyr::desc(in_target_range), distance_to_target, dplyr::desc(threshold))

  chosen_threshold <- summary$threshold[1]
  idx <- which(vapply(trials, function(x) isTRUE(all.equal(x$threshold, chosen_threshold)), logical(1)))[1]
  chosen <- trials[[idx]]
  chosen$trial_summary <- summary
  chosen
}
