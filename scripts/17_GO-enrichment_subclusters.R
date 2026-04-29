#GO-term analysis from subcluster-based fisher-consensus results
#after this analysis, GO-term results were summarized with revigo (revigo.irb.hr)

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(dplyr)
library(rrvgo)
library(cowplot)
library(patchwork)

base_output <- file.path("~/HAECA/global_analysis/fisher-consensus_analysis_subclusters")
setwd(base_output)
GO_output <- file.path(paste0(base_output, "/GO_analysis"))
dir.create(GO_output, showWarnings = T)

#define clusters to analyse
clusters <- c("artery","capillary","vein")

run_GO_analysis <- function(GO_genes, dir, lfc_cutoff = 0, top_n = NULL, 
                            pval_cutoff = 0.05, pval_max = 0.1, 
                            min_genes = 10) {
  
  if (dir == "up") {
    genes_dir <- GO_genes %>% filter(log2FoldChange > 0)
  } else {
    genes_dir <- GO_genes %>% filter(log2FoldChange < 0)
  }
  
  take_top_n <- function(df, n, col_to_pull = "gene") {
    n_available <- nrow(df)
    n_used <- min(n_available, n)
    genes <- df %>% head(n = n_used) %>% pull(!!sym(col_to_pull))
    return(list(genes = genes, n_used = n_used, n_available = n_available))
  }
  
  valid_data <- genes_dir %>%
      filter(!is.na(pvalue))
    
    sig_data <- valid_data %>%
      filter(pvalue < pval_max) %>%
      arrange(pvalue)
    
    nonsig_data <- valid_data %>%
      filter(pvalue >= pval_max) %>%
      arrange(desc(abs(log2FoldChange)))
    
    combined <- bind_rows(sig_data, nonsig_data)
    
    result <- take_top_n(combined, top_n)
    
    n_sig_used <- min(nrow(sig_data), result$n_used)
    n_nonsig_used <- max(0, result$n_used - nrow(sig_data))
    
    GO_genes_filtered <- result$genes
    suffix <- paste0("top", top_n, "_bypval")
    
    description <- paste0(
      "Top ", result$n_used, " (",
      n_sig_used, " p<", pval_max, " + ",
      n_nonsig_used, " ≥p)",
      ifelse(result$n_used < top_n,
             paste0(" [", result$n_available, " avail]"), "")
    )
  
  return(list(
    genes = GO_genes_filtered,
    suffix = suffix,
    description = description,
    n_genes = length(GO_genes_filtered),
    n_requested = ifelse(is.null(top_n), NA, top_n)
  ))
}
  
#set correct directory
method_dir <- file.path(paste0(base_output, "/consensus_median_subclusters"))

plots_list <- list()
  
  for(cluster in clusters) {
    
    #set directory for cluster
    cluster_dir <- paste0(method_dir, "/", cluster)
    
    # Check if combined results file exists
    combined_file <- paste0(cluster_dir, "/median_Combined_results.csv")
    if (!file.exists(combined_file)) {
      message("Skipping ", cluster, ": Combined results file not found")
      next
    }
    
    #load data
    final_combined <- read.csv(combined_file)
    
    go_file <- paste0(base_output, "/plots_median/genes_consensus_shared_median_", cluster, ".csv")
    
    if (!file.exists(go_file)) {
      message("Skipping ", cluster, ": file not found")
      # Add empty plots for this cluster
      plot_list[[paste0(cluster, "_up")]] <- NULL
      plot_list[[paste0(cluster, "_down")]] <- NULL
      next
    }
    
    GO_genes <- read.csv(go_file)
    
    # Background genes
    genes_all <- final_combined %>%
      distinct(gene) %>%
      pull()
    
    ref <- bitr(genes_all, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db) %>%
      distinct(SYMBOL, .keep_all = TRUE) %>%
      filter(SYMBOL %in% genes_all)
    
  for (dir in c('up', 'down')) {
    
    # Get filtered genes
    result <- run_GO_analysis(
      GO_genes = GO_genes,
      dir = dir,
      top_n = top_n,
      pval_cutoff = 0.05,
      pval_max = pval_max_threshold
    )
    
    GO_genes_filtered <- result$genes
    suffix <- result$suffix
    description <- result$description
    n_input_genes <- result$n_genes
    
    strategy_suffix <- suffix
    strategy_description <- description

        # Skip if no genes
    if (length(GO_genes_filtered) == 0) {
      message("No genes for ", dir, " ", suffix, ". Skipping.")
      next
    }
    
    # Convert gene symbols
    sig_gene_ids <- bitr(GO_genes_filtered, fromType = "SYMBOL",
                         toType = c("ENTREZID", "ENSEMBL"),
                         OrgDb = org.Hs.eg.db) %>%
      distinct(SYMBOL, .keep_all = TRUE)
    
    # Run GO enrichment
    GOres <- enrichGO(
      gene = sig_gene_ids$SYMBOL,
      universe = ref$SYMBOL,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pvalueCutoff = 0.25,
      qvalueCutoff = 0.25
    )
    
    # Save GO results
    write.csv(as.data.frame(GOres),
              paste0(GO_output, "/GO_result_", cluster, "_", dir, "_", strategy_suffix, ".csv"),
              row.names = F)
    
    # Skip if no results
    if (is.null(GOres) || nrow(as.data.frame(GOres)) == 0) {
      message("No GO results for ", dir, " ", suffix, ".")
      next
    }
    
  }
    
  }
