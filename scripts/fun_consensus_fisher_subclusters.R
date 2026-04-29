
read_genestats <- function(stats_dir,
                           pattern = "^DESeq2_gene_stability_.*_subsampling\\.csv$") {
  files <- list.files(stats_dir, pattern = pattern, full.names = TRUE)
  map_df(files, function(f) {
    df <- read_csv(f, show_col_types = FALSE)
    df
  })
}

read_real_deseq <- function(real_dir,
                            pattern = "^DESeq2_.*_real\\.csv$") {
  files <- list.files(real_dir, pattern = pattern, full.names = TRUE)
  
  map_df(files, function(f) {
    df <- read_csv(f, show_col_types = FALSE)
    
    fname <- basename(f)
    core <- fname %>%
      str_remove("^DESeq2_") %>%
      str_remove("_real\\.csv$")
    
    parts <- str_split(core, "_", simplify = TRUE)
    tissue_from_file <- parts[1]
    cluster_from_file <- paste(parts[-1], collapse = "_")
    
    df %>%
      mutate(
        # Use existing columns if present, otherwise extract from filename
        tissue = if ("tissue" %in% names(df)) tissue else tissue_from_file,
        cluster = if ("cluster" %in% names(df)) cluster else cluster_from_file
      )
  }) %>%
    filter(!is.na(pvalue), !is.na(log2FoldChange)) %>%
    mutate(
      direction = case_when(
        log2FoldChange > 0 ~ "up",
        log2FoldChange < 0 ~ "down",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(direction))
}

compute_consensus_deg <- function(all_stats,
                                  prop_sig_cutoff = 0.60,
                                  min_tissues = 3,
                                  min_abs_median_lfc = 0.25) {
  
  robust <- all_stats %>%
    filter(prop_sig >= prop_sig_cutoff) %>%
    mutate(
      direction = case_when(
        mean_lfc > 0 ~ "up",
        mean_lfc < 0 ~ "down",
        TRUE ~ NA_character_
      ),
      stability_index = 1 / (1 + sd_lfc)
    ) %>%
    filter(!is.na(direction))
  
  consensus <- robust %>%
    group_by(cluster, gene, direction) %>%
    summarise(
      n_tissues = n_distinct(tissue),
      mean_effect = mean(mean_lfc, na.rm = TRUE),
      median_effect = median(mean_lfc, na.rm = TRUE),
      mean_stability = mean(stability_index, na.rm = TRUE),
      tissues = paste(sort(unique(tissue)), collapse = "; "),
      .groups = "drop"
    )
  
  consensus %>%
    filter(
      n_tissues >= min_tissues,
      abs(median_effect) >= min_abs_median_lfc
    ) %>%
    transmute(
      cluster,
      gene,
      consensus_direction = direction,
      consensus_n_tissues = n_tissues,
      consensus_mean_lfc = mean_effect,
      consensus_median_lfc = median_effect,
      consensus_tissues = tissues,
      mean_stability
    )
}

run_aw_fisher <- function(real_all, min_tissues = 2, cluster_name = NULL) {
  
  if (is.null(cluster_name)) {
    cluster_name <- unique(real_all$cluster)[1]
    if (is.na(cluster_name) || is.null(cluster_name)) {
      cluster_name <- "unknown"
    }
  }
  
  n_tiss <- dplyr::n_distinct(real_all$tissue)
  
  if (n_tiss < min_tissues) {
    return(list(
      ok = FALSE,
      reason = paste0("Only ", n_tiss, " tissue(s) available; need >= ", min_tissues),
      aw_obj = NULL
    ))
  }
  
  gene_tissue_counts <- real_all %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(n_tissues = dplyr::n_distinct(tissue), .groups = "drop")
  
  genes_sufficient <- gene_tissue_counts %>%
    dplyr::filter(n_tissues >= min_tissues) %>%
    dplyr::pull(gene)
  
  if (length(genes_sufficient) == 0) {
    return(list(
      ok = FALSE,
      reason = paste0("No genes present in >= ", min_tissues, " tissues"),
      aw_obj = NULL
    ))
  }
  
  real_all_aw <- real_all %>%
    dplyr::filter(gene %in% genes_sufficient) %>%
    dplyr::arrange(gene, tissue)
  
  if (nrow(real_all_aw) == 0) {
    return(list(ok = FALSE, reason = "No genes after filtering", aw_obj = NULL))
  }
  
  aw_res <- tryCatch({
    real_all_aw %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(
        tissues = list(tissue),
        p_matrix = list(matrix(pvalue, nrow = 1)),
        aw = list(AWFisher::AWFisher_pvalue(p_matrix[[1]])),
        p_aw = aw[[1]]$pvalues,
        weights = list(aw[[1]]$weights),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        cluster = cluster_name,
        p_aw_fdr = stats::p.adjust(p_aw, method = "BH")
      )
  }, error = function(e) {
    message("AWFisher error: ", e$message)
    return(NULL)
  })
  
  if (is.null(aw_res)) {
    return(list(ok = FALSE, reason = "AWFisher computation failed", aw_obj = NULL))
  }

  weights_long <- aw_res %>%
    dplyr::select(cluster, gene, tissues, weights, p_aw, p_aw_fdr) %>%
    dplyr::mutate(
      tissues = purrr::map(tissues, as.character),
      weights = purrr::map(weights, as.numeric)
    ) %>%
    tidyr::unnest(c(tissues, weights)) %>%
    dplyr::rename(tissue = tissues)
  
  consistency_table <- weights_long %>%
    dplyr::group_by(cluster, gene) %>%
    dplyr::summarise(
      n_tissues = dplyr::n(),
      n_selected = sum(weights == 1),
      prop_selected = n_selected / n_tissues,
      pattern = dplyr::case_when(
        n_selected == n_tissues ~ "globally_consistent",
        n_selected == 1 ~ "tissue_specific",
        TRUE ~ "subset_consistent"
      ),
      p_aw = dplyr::first(p_aw),
      p_aw_fdr = dplyr::first(p_aw_fdr),
      .groups = "drop"
    )
  
  direction_consistency <- real_all_aw %>%
    dplyr::mutate(cluster = cluster_name) %>%
    dplyr::group_by(cluster, gene) %>%
    dplyr::summarise(
      n_up = sum(log2FoldChange > 0),
      n_down = sum(log2FoldChange < 0),
      direction_consistent = (n_up == 0 | n_down == 0),
      mean_lfc = mean(log2FoldChange, na.rm = TRUE),
      .groups = "drop"
    )
  
  aw_meta <- consistency_table %>%
    dplyr::left_join(direction_consistency, by = c("cluster", "gene")) %>%
    dplyr::mutate(
      aw_direction = dplyr::if_else(mean_lfc > 0, "up", "down"),
      aw_direction_consistent = direction_consistent
    ) %>%
    dplyr::select(
      cluster, gene, p_aw, p_aw_fdr, pattern,
      aw_n_tissues_total = n_tissues,
      aw_n_tissues_selected = n_selected,
      prop_selected,
      aw_direction, aw_direction_consistent, mean_lfc
    )
  
  selected_tissues_df <- weights_long %>%
    dplyr::filter(weights == 1) %>%
    dplyr::group_by(cluster, gene) %>%
    dplyr::summarise(
      aw_selected_tissues = paste(tissue, collapse = "; "),
      .groups = "drop"
    )
  
  list(
    ok = TRUE,
    reason = NA_character_,
    aw_obj = list(
      aw_res = aw_res,
      weights_long = weights_long,
      consistency_table = consistency_table,
      direction_consistency = direction_consistency,
      aw_meta = aw_meta,
      selected_tissues = selected_tissues_df
    )
  )
}

combine_consensus_and_aw <- function(consensus_deg, aw_obj) {
  
  consensus_deg %>%
    full_join(aw_obj$aw_meta, by = c("cluster", "gene")) %>%
    left_join(aw_obj$selected_tissues, by = c("cluster", "gene")) %>%
    mutate(
      direction_agreement = case_when(
        is.na(consensus_direction) | is.na(aw_direction) ~ "insufficient_data",
        consensus_direction == aw_direction ~ "agree",
        TRUE ~ "conflict"
      ),
      final_direction = if_else(direction_agreement == "agree",
                                consensus_direction, NA_character_)
    ) %>%
    arrange(cluster, p_aw_fdr)
}

run_fisher_pipeline <- function(stats_dir, real_dir, out_base, prefix,
                                prop_sig_cutoff = 0.40,
                                min_tissues_consensus = 2,
                                min_abs_median_lfc = 0.2,
                                min_tissues_aw = 2) {
  
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
  
  all_stats <- read_genestats(stats_dir)
  real_all  <- read_real_deseq(real_dir)
  
  if (nrow(all_stats) == 0 || nrow(real_all) == 0) {
    stop("Input data is empty. Check stats_dir and real_dir paths.")
  }
  
  clusters <- sort(unique(intersect(all_stats$cluster, real_all$cluster)))
  
  if (length(clusters) == 0) {
    warning("No overlapping clusters between stats and real data.")
    return(list())
  }
  
  res <- setNames(vector("list", length(clusters)), clusters)
  
  for (cl in clusters) {
    message("\n=== META for cluster: ", cl, " ===")
    
    stats_cl <- all_stats %>% filter(cluster == cl)
    real_cl  <- real_all  %>% filter(cluster == cl)
    
    tissues_ok <- intersect(unique(stats_cl$tissue), unique(real_cl$tissue))
    
    if (length(tissues_ok) < min_tissues_aw) {
      message("Skipping ", cl, ": only ", length(tissues_ok),
              " tissue(s) available; need >= ", min_tissues_aw, " for AW-Fisher.")
      res[[cl]] <- list(
        cluster = cl,
        skipped = TRUE,
        reason = paste0("Only ", length(tissues_ok), " tissues available")
      )
      next
    }
    
    stats_cl <- stats_cl %>% filter(tissue %in% tissues_ok)
    real_cl  <- real_cl  %>% filter(tissue %in% tissues_ok)
  
    if (nrow(real_cl) == 0) {
      message("Skipping ", cl, ": No genes remaining after tissue filtering.")
      next
    }
    
    cl_safe <- gsub("[^A-Za-z0-9._-]+", "_", cl)
    out_dir <- file.path(out_base, cl_safe)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    consensus <- compute_consensus_deg(
      stats_cl,
      prop_sig_cutoff = prop_sig_cutoff,
      min_tissues = min_tissues_consensus,
      min_abs_median_lfc = min_abs_median_lfc
    )
    
    aw_out <- run_aw_fisher(real_cl, min_tissues = min_tissues_aw, cluster_name = cl)
    
    if (!aw_out$ok) {
      message("Skipping ", cl, ": ", aw_out$reason)
      res[[cl]] <- list(
        cluster = cl,
        skipped = TRUE,
        reason = aw_out$reason,
        consensus_deg = consensus
      )
      next
    }
    
    combined <- combine_consensus_and_aw(consensus, aw_out$aw_obj)
    
    readr::write_csv(combined, file.path(out_dir, paste0(prefix, "_Combined_results.csv")))
    readr::write_csv(
      aw_out$aw_obj$consistency_table %>% 
        left_join(aw_out$aw_obj$direction_consistency, by = "gene"),
      file.path(out_dir, paste0(prefix, "_AWFisher_results.csv"))
    )
    
    res[[cl]] <- list(
      cluster = cl,
      skipped = FALSE,
      tissues_used = tissues_ok,
      consensus_deg = consensus,
      aw = aw_out$aw_obj,
      final_combined = combined
    )
  }
  
  invisible(res)
}

compare_results <- function(res_median_by_cluster,
                                       res_cutoff_by_cluster,
                                       out_base,
                                       name_a = "median",
                                       name_b = "cutoff") {
  
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
  
  safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
  
  clusters <- sort(unique(c(names(res_median_by_cluster), names(res_cutoff_by_cluster))))
  
  overall_summary <- list()
  
  for (cl in clusters) {
    a <- res_median_by_cluster[[cl]]$final_combined
    b <- res_cutoff_by_cluster[[cl]]$final_combined
    if (is.null(a) && is.null(b)) next
    
    cl_safe <- safe_name(cl)
    out_dir <- file.path(out_base, cl_safe)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    genes_a <- if (!is.null(a)) {
      a %>% filter(!is.na(final_direction)) %>%
        select(gene, final_direction) %>% distinct() %>%
        rename(!!paste0("final_direction_", name_a) := final_direction)
    } else {
      tibble(gene = character(),
             !!paste0("final_direction_", name_a) := character())
    }
    
    genes_b <- if (!is.null(b)) {
      b %>% filter(!is.na(final_direction)) %>%
        select(gene, final_direction) %>% distinct() %>%
        rename(!!paste0("final_direction_", name_b) := final_direction)
    } else {
      tibble(gene = character(),
             !!paste0("final_direction_", name_b) := character())
    }
    
    comp <- full_join(genes_a, genes_b, by = "gene") %>%
      mutate(
        cluster = cl,
        in_a = !is.na(.data[[paste0("final_direction_", name_a)]]),
        in_b = !is.na(.data[[paste0("final_direction_", name_b)]]),
        shared = in_a & in_b,
        status = case_when(
          shared & (.data[[paste0("final_direction_", name_a)]] ==
                      .data[[paste0("final_direction_", name_b)]]) ~ "shared_same_direction",
          shared ~ "shared_different_direction",
          in_a ~ paste0(name_a, "_only"),
          in_b ~ paste0(name_b, "_only"),
          TRUE ~ "unknown"
        )
      ) %>%
      arrange(desc(status == "shared_same_direction"), gene)
    
    summ <- comp %>%
      count(cluster, status, name = "n") %>%
      group_by(cluster) %>%
      mutate(pct = round(100 * n / sum(n), 1)) %>%
      ungroup()
    
    write_csv(comp, file.path(out_dir, paste0(name_a, "_vs_", name_b, "_", cl_safe, ".csv")))
    write_csv(summ, file.path(out_dir, paste0("summary", name_a, "_vs_", name_b, "_", cl_safe, ".csv")))
    
    overall_summary[[cl]] <- summ
  }
  
  overall_summary_df <- bind_rows(overall_summary)
  write_csv(overall_summary_df, file.path(out_base, paste0("summary_all_clusters_", name_a, "_vs_", name_b, ".csv")))
  
  invisible(overall_summary_df)
}