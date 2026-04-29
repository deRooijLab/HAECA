#### Load libraries ####
library(tidyverse)
library(readr)
library(AWFisher)
library(ComplexHeatmap)
library(circlize)
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)
library(tidyplots)

# Directories
base_input <- file.path("~/HAECA/global_analysis/downsampling")
base_output <- ("~/HAECA/global_analysis/consensus_analysis")
dir.create(base_output, showWarnings = TRUE)
plot_output <- (paste0(base_output, "/plots"))
dir.create(plot_output, showWarnings = TRUE)

source(paste0(base_output, '/fun_consensus_fisher.R'))

# Run median
res_median <- run_fisher_pipeline(
  
  stats_dir = file.path(base_input, "output_median_parallel"),
  real_dir = file.path(base_input, "output_median_parallel"),
  output_dir = file.path(base_output, "output_median"),
  prefix = "median"
)

# Run cutoff
res_cutoff <- run_fisher_pipeline(
  stats_dir = file.path(base_input, "output_cutoff_parallel"),
  real_dir = file.path(base_input, "output_cutoff_parallel"),
  output_dir = file.path(base_output, "output_cutoff"),
  prefix = "cutoff"
)

########## meta-analsis #########
color_gradient <- c("#D1422F","#F4AB5C","#e9c46a","#ACC8BE","#1A5B5B","#264653")

#prepare function
plot_consensus <- function(top_gene_data, method = "", title = paste0("Top up- and down-regulated genes (consensus) for ", method)) {
  
  #save top genes with separate UP and DOWN columns
  genes_up <- top_gene_data %>% filter(log2FoldChange > 0) %>% pull(gene)
  genes_dn <- top_gene_data %>% filter(log2FoldChange < 0) %>% pull(gene)
  
  #pad to equal length for data frame
  max_len <- max(length(genes_up), length(genes_dn))
  genes_df <- data.frame(
    UP   = c(genes_up, rep(NA, max_len - length(genes_up))),
    DOWN = c(genes_dn, rep(NA, max_len - length(genes_dn)))
  )
  write.csv(genes_df, paste0("top_genes_consensus_up_dn_", method, ".csv"),
            quote = FALSE, row.names = FALSE, na = "")
  
  #combined dotplot (both up & down) as lollipop
  plot_df <- top_gene_data %>%
    mutate(
      p_safe = pmax(pvalue, .Machine$double.xmin),
      neglog10p = -log10(p_safe),
      gene_ord = fct_reorder(gene, log2FoldChange)
    )
  
  #spacing for labels
  rng <- range(plot_df$log2FoldChange, na.rm = TRUE)
  pad <- 0.08 * diff(rng)
  
  plot_df <- plot_df %>%
    mutate(
      x_lab = if_else(log2FoldChange >= 0, log2FoldChange + pad, log2FoldChange - pad),
      hjust_lab = if_else(log2FoldChange >= 0, 0, 1)
    )
  
  p <- ggplot(plot_df, aes(x = log2FoldChange, y = gene_ord)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey70", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = log2FoldChange, yend = gene_ord),
                 color = "grey30", linewidth = 0.7) +
    geom_point(aes(fill = neglog10p, size = consensus_n_tissues),
               shape = 21, color = "black", stroke = 0.25, alpha = 1) +
    geom_text(aes(x = x_lab, label = gene, hjust = hjust_lab),
              size = 3.2) +
    scale_y_discrete(labels = NULL) +
    scale_fill_gradientn(colors = rev(color_gradient), name = expression(-log[10](p))) +
    scale_size_continuous(
      range = c(6, 12),
      breaks = sort(unique(plot_df$consensus_n_tissues)),
      name = "# tissues"
    ) +
    scale_x_continuous(
      limits = c(rng[1] - 3*pad, rng[2] + 3*pad),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0("Top up- and down-regulated genes (consensus) for ", method),
      x = "Median consensus log2FC",
      y = "Gene"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.35),
      panel.grid.minor.x = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(5.5, 40, 5.5, 40)
    )
  
  return(p)
  
}

#take dataframes from consensus/fisher test for plotting
final_median <- res_median$final_combined
final_cutoff <- res_cutoff$final_combined

final_list <- list(
  median = final_median,
  cutoff = final_cutoff
)

for(i in 1:length(final_list)) {
  
  method <- names(final_list)[i]
  final_combined <- final_list[[i]]

####plot top genes across tissues####
#pick values based on consensus analysis (coming from permutations)
#exclude tissue-specific genes, select shared direction in fisher and consensus

GO_genes <- final_combined %>%
  filter(
    direction_agreement == "agree" &
      pattern != 'tissue_specific' &
      !is.na(consensus_median_lfc))  %>%
  mutate(log2FoldChange = consensus_median_lfc,
         pvalue = p_aw_fdr)

write.csv(GO_genes, paste0(base_output, "/genes_consensus_shared_", method, ".csv"), row.names = F)

#selecting genes based on lowest pvalue
top_gene_data <- GO_genes %>%
  dplyr::select(gene, pvalue, log2FoldChange, consensus_n_tissues) %>%
  arrange(log2FoldChange) %>%
  { bind_rows(
    filter(., log2FoldChange < 0) %>% arrange(pvalue) %>% slice_head(n = 25),
    filter(., log2FoldChange > 0) %>% arrange(pvalue) %>% slice_head(n = 25)
  ) } %>%
  distinct(gene, .keep_all = TRUE)

p <- plot_consensus(top_gene_data, method = method)
ggsave(paste0(plot_output, "/top25_consensus_sig_both_", method, ".pdf"), plot = p,
       width = 6, height = 12)

#export data for supp. table S6 and S8
write.csv(top_gene_data, paste0(plot_output, "/top25_consensus_sig_both_", method, ".csv"))

}