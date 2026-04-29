#### Load libraries ####
library(tidyverse)
library(readr)
library(AWFisher)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(viridis)
library(ggplot2)
library(pheatmap)

# Directories
base_input <- file.pathfile.path("~/HAECA/global_analysis/downsampling_subclusters")
base_output <- file.path("~/HAECA/global_analysis/consensus_analysis_subclusters")
setwd(base_output)

source("fun_consensus_fisher_subclusters.R")

res_median <- run_fisher_pipeline(
  stats_dir = file.path(base_input, "output_median"),
  real_dir  = file.path(base_input, "output_median"),
  out_base  = file.path(base_output, "consensus_median_subclusters"),
  prefix    = "median"
)

make_heatmap <- function(heatmap_data, title){
  
  heatmap_input <- heatmap_data %>%
    dplyr::select(consensus_median_lfc, cluster, gene) %>%
    pivot_wider(
      names_from = cluster,
      values_from = consensus_median_lfc,
      values_fill = NA
    ) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  # heatmap_matrix_log <- -log10(heatmap_input)
  heatmap_matrix_log <- heatmap_input
  heatmap_matrix_log[is.infinite(heatmap_matrix_log)] <- 10  # Cap Inf
  heatmap_matrix_log[is.na(heatmap_matrix_log)] <- 0
  
  color_hm <- if (dir == 'up') viridis(100) else rev(viridis(100))
  
  hm <- pheatmap(t(heatmap_matrix_log),
                 main = title,
                 color = color_hm,
                 border_color = "NA",
                 na_col = "grey90",
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 fontsize_row = 7,
                 fontsize_col = 10,
                 cellwidth = 10,
                 cellheight = 10)
  
  return(hm)
}

#plotting heatmaps for vascular beds
#set directory
method_dir <- file.path(paste0(base_output, "/fisher_median_cluster"))
output_dir <- file.path(paste0(base_output, "/heatmaps_", method))
dir.create(output_dir)
  
#pick available clusters
clusters <- list.dirs(method_dir, full.names = FALSE, recursive = FALSE)
  
# Accumulate data from all clusters into one dataframe
data_all <- map_df(clusters, function(cluster) {
    message("Processing: ", method, " - ", cluster)
    
    cluster_dir <- paste0(method_dir, "/", cluster)
    combined_file <- paste0(cluster_dir, "/", method, "_Combined_results.csv")
    
    if (!file.exists(combined_file)) {
      message("Skipping ", cluster, ": Combined results file not found")
      return(NULL)
    }
    
    df <- read.csv(combined_file)
    
    if (nrow(df) == 0) {
      message("Skipping ", cluster, ": No data")
      return(NULL)
    }
    
    df <- df %>%
      mutate(cluster = cluster)
    
    return(df)
  })
  
#filtering for selection of targets for heatmap
  top_gene_data <- data_all %>%
    filter(cluster %in% c('capillary', 'artery', 'vein') &
             direction_agreement == 'agree' &
             pattern != 'tissue_specific')
  
  for (dir in c("up", "down")) {
    
    heatmap_data <- top_gene_data %>%
      dplyr::select(gene, cluster, p_aw_fdr, consensus_median_lfc, consensus_n_tissues) %>%
      filter(if (dir == "up") consensus_median_lfc > 0 else consensus_median_lfc < 0) %>%
      arrange(p_aw_fdr) %>%
      head(n = 80)

    #Figure 2f+j: group genes by cluster, then select top 20 genes each that are shared in at least 3 tissues
    heatmap_data <- top_gene_data %>%
      dplyr::select(gene, cluster, p_aw_fdr, consensus_median_lfc, consensus_n_tissues) %>%
      group_by(cluster) %>%
      filter(if (dir == "up") consensus_median_lfc > 0 else consensus_median_lfc < 0) %>%
      filter(consensus_n_tissues >= 3) %>%
      arrange(p_aw_fdr) %>%
      slice_head(n = 20) %>%
      ungroup()
    
    hm1 <- make_heatmap(heatmap_data,
                        title = paste0(dir, "-regulated genes \n AWFisher consensus median lfc across clusters, \n
                            top 20 genes per cluster, shared in >= 3 tissues, ranked by pvalue"))
    
    #group by cluster, then select top 35 genes each that are shared in at least 3 tissues
    heatmap_data <- top_gene_data %>%
      dplyr::select(gene, cluster, p_aw_fdr, consensus_median_lfc, consensus_n_tissues) %>%
      group_by(cluster) %>%
      filter(if (dir == "up") consensus_median_lfc > 0 else consensus_median_lfc < 0) %>%
      filter(consensus_n_tissues >= 3) %>%
      arrange(p_aw_fdr) %>%
      slice_head(n = 35) %>%
      ungroup()
    
    hm2 <- make_heatmap(heatmap_data,
                        title = paste0(dir, "-regulated genes \n AWFisher consensus median lfc across clusters, \n
                            top 35 genes per cluster, shared in >= 3 tissues, ranked by pvalue"))
    
    
    #save all
    pdf(file = paste0(output_dir, "/heatmaps_", dir, ".pdf"), width = 18, height = 7)
    for(n in c(1:2)){
      hm <- get((paste0("hm",n)))
      print(hm)
    }
    dev.off()
    
  }