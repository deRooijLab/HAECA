
read_genestats <- function(stats_dir, pattern = "^GeneStats_.*_allRepeats.csv$") {
  list.files(stats_dir, pattern = pattern, full.names = TRUE) %>%
    map_df(function(f) {
      tissue <- basename(f) %>% 
        str_remove("^GeneStats_") %>% 
        str_remove("_allRepeats.csv$")
      read_csv(f, show_col_types = FALSE) %>% 
        mutate(tissue = tissue)
    })
}

read_real_deseq <- function(real_dir, pattern = "^DESeq2_.*_real.csv$") {
  list.files(real_dir, pattern = pattern, full.names = TRUE) %>%
    map_df(function(f) {
      tissue <- basename(f) %>% 
        str_remove("^DESeq2_") %>% 
        str_remove("_real.csv$")
      read_csv(f, show_col_types = FALSE) %>% 
        mutate(tissue = tissue)
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
  robust_by_tissue <- all_stats %>%
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
  
  consensus <- robust_by_tissue %>%
    group_by(gene, direction) %>%
    summarise(
      n_tissues = n_distinct(tissue),
      mean_effect = mean(mean_lfc),
      median_effect = median(mean_lfc),
      mean_stability = mean(stability_index),
      tissues = paste(sort(unique(tissue)), collapse = "; "),
      .groups = "drop"
    )

  consensus %>%
    filter(
      n_tissues >= min_tissues,
      abs(median_effect) >= min_abs_median_lfc
    ) %>%
    transmute(
      gene,
      consensus_direction = direction,
      consensus_n_tissues = n_tissues,
      consensus_mean_lfc = mean_effect,
      consensus_median_lfc = median_effect,
      consensus_tissues = tissues,
      mean_stability
    )
}

run_aw_fisher <- function(real_all, min_tissues = 2) {
  real_all_aw <- real_all %>%
    group_by(gene) %>%
    filter(n() >= min_tissues) %>%
    ungroup() %>%
    arrange(gene, tissue)
  
  aw_res <- real_all_aw %>%
    group_by(gene) %>%
    summarise(
      tissues = list(tissue),
      p_matrix = list(matrix(pvalue, nrow = 1)),
      aw = list(AWFisher_pvalue(p_matrix[[1]])),
      p_aw = aw[[1]]$pvalues,
      weights = list(aw[[1]]$weights),
      .groups = "drop"
    ) %>%
    mutate(p_aw_fdr = p.adjust(p_aw, method = "BH"))
  
  weights_long <- aw_res %>%
    select(gene, tissues, weights, p_aw, p_aw_fdr) %>%
    mutate(
      tissues = map(tissues, as.character),
      weights = map(weights, as.numeric)
    ) %>%
    unnest(c(tissues, weights))
  
  consistency_table <- weights_long %>%
    group_by(gene) %>%
    summarise(
      n_tissues = n(),
      n_selected = sum(weights == 1),
      prop_selected = n_selected / n_tissues,
      pattern = case_when(
        n_selected == n_tissues ~ "globally_consistent",
        n_selected == 1 ~ "tissue_specific",
        TRUE ~ "subset_consistent"
      ),
      p_aw = first(p_aw),
      p_aw_fdr = first(p_aw_fdr),
      .groups = "drop"
    )
  
  direction_consistency <- real_all %>%
    group_by(gene) %>%
    summarise(
      n_up = sum(log2FoldChange > 0),
      n_down = sum(log2FoldChange < 0),
      direction_consistent = (n_up == 0 | n_down == 0),
      mean_lfc = mean(log2FoldChange),
      .groups = "drop"
    )
  
  aw_meta <- consistency_table %>%
    left_join(direction_consistency, by = "gene") %>%
    mutate(
      aw_direction = if_else(mean_lfc > 0, "up", "down"),
      aw_direction_consistent = direction_consistent
    ) %>%
    select(
      gene, p_aw, p_aw_fdr, pattern,
      aw_n_tissues_total = n_tissues,
      aw_n_tissues_selected = n_selected,
      prop_selected,
      aw_direction, aw_direction_consistent, mean_lfc
    )
  
  selected_tissues_df <- weights_long %>%
    filter(weights == 1) %>%
    group_by(gene) %>%
    summarise(aw_selected_tissues = paste(tissues, collapse = "; "), .groups = "drop")
  
  list(
    aw_res = aw_res,
    weights_long = weights_long,
    consistency_table = consistency_table,
    direction_consistency = direction_consistency,
    aw_meta = aw_meta,
    selected_tissues = selected_tissues_df
  )
}

combine_consensus_and_aw <- function(consensus_deg, aw_obj) {
  consensus_deg %>%
    full_join(aw_obj$aw_meta, by = "gene") %>%
    left_join(aw_obj$selected_tissues, by = "gene") %>%
    mutate(
      direction_agreement = case_when(
        is.na(consensus_direction) | is.na(aw_direction) ~ "insufficient_data",
        consensus_direction == aw_direction ~ "agree",
        TRUE ~ "conflict"
      ),
      final_direction = if_else(direction_agreement == "agree", 
                                consensus_direction, NA_character_)
    ) %>%
    arrange(p_aw_fdr)
}

run_fisher_pipeline <- function(stats_dir, 
                                real_dir, 
                                output_dir, 
                                prefix,
                                prop_sig_cutoff = 0.60,
                                min_tissues_consensus = 3,
                                min_abs_median_lfc = 0.25,
                                min_tissues_aw = 2) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 1. Read data
  
  all_stats <- read_genestats(stats_dir)
  real_all <- read_real_deseq(real_dir)
  
  
  # 2. Consensus DEGs
  
  consensus_deg <- compute_consensus_deg(
    all_stats,
    prop_sig_cutoff = prop_sig_cutoff,
    min_tissues = min_tissues_consensus,
    min_abs_median_lfc = min_abs_median_lfc
  )
  
  # 3. AW-Fisher
  
  aw_obj <- run_aw_fisher(real_all, min_tissues = min_tissues_aw)
  
  # 4. Combine
  
  final_combined <- combine_consensus_and_aw(consensus_deg, aw_obj)
  
  # 5. Save
  
  write_csv(
    aw_obj$consistency_table %>% left_join(aw_obj$direction_consistency, by = "gene"),
    file.path(output_dir, paste0(prefix, "_AWFisher_results.csv"))
  )
  write_csv(
    final_combined,
    file.path(output_dir, paste0(prefix, "_Combined_results.csv"))
  )
  
  # Return everything
  
  list(
    all_stats = all_stats,
    real_all = real_all,
    consensus_deg = consensus_deg,
    aw = aw_obj,
    final_combined = final_combined
  )
}