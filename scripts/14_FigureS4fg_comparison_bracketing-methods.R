#comparison of genes coming from median-based and cutoff-based age-bracket methods
#Figure S4f & S4g

library(tidyverse)
library(readr)
library(ggrepel)
library(viridis)

base_input <- file.path("~/HAECA/global_analysis/consensus_analysis")
base_output <- file.path("~/HAECA/global_analysis/consensus_comparison")
setwd(base_output)

output_dir <- file.path(paste0(base_output, "/plots"))
dir.create(output_dir)

# Load combined data from both methods
data_combined <- map_df(c("median", "cutoff"), function(method) {
  method_dir <- file.path(paste0(base_input, "/output_", method))
  file_path <- paste0(method_dir, "/", method, "_Combined_results.csv")
  if (!file.exists(file_path)) {
    message("File not found: ", file_path)
    return(NULL)
  }
  read.csv(file_path) %>%
    mutate(method = method)
})

#combine and pre-select for genes represented in >= 2 tissues
comparison_data_both <- data_combined %>%
  filter(pattern != 'tissue_specific' & !is.na(final_direction)) %>%
  dplyr::select(gene, method, final_direction, consensus_median_lfc, p_aw_fdr, mean_lfc,
                consensus_n_tissues, aw_n_tissues_total) %>%
  pivot_wider(
    id_cols = c(gene, final_direction),
    names_from = method,
    values_from = c(consensus_median_lfc, p_aw_fdr, mean_lfc,
                    consensus_n_tissues, aw_n_tissues_total),
    names_glue = "{.value}_{method}"
  ) %>%
  mutate(
    # Categorize genes
    in_both = !is.na(consensus_median_lfc_median) & !is.na(consensus_median_lfc_cutoff),
    only_median = !is.na(consensus_median_lfc_median) & is.na(consensus_median_lfc_cutoff),
    only_cutoff = is.na(consensus_median_lfc_median) & !is.na(consensus_median_lfc_cutoff),
    # Category as factor for plotting
    category = case_when(
      in_both ~ "Both methods",
      only_median ~ "Median only",
      only_cutoff ~ "Cutoff only"
    ),
    # P-value summaries
    min_pval = pmin(p_aw_fdr_median, p_aw_fdr_cutoff, na.rm = TRUE),
    max_pval = pmax(p_aw_fdr_median, p_aw_fdr_cutoff, na.rm = TRUE),
    # Significance flags
    sig_median = !is.na(p_aw_fdr_median) & p_aw_fdr_median < 0.05,
    sig_cutoff = !is.na(p_aw_fdr_cutoff) & p_aw_fdr_cutoff < 0.05,
    sig_both = in_both & sig_median & sig_cutoff,
    # LFC comparisons
    lfc_diff = consensus_median_lfc_median - consensus_median_lfc_cutoff,
    lfc_mean = (consensus_median_lfc_median + consensus_median_lfc_cutoff) / 2,
    # Direction consistency between methods
    consistent_lfc = case_when(
      !in_both ~ NA,
      consensus_median_lfc_median > 0 & consensus_median_lfc_cutoff > 0 ~ TRUE,
      consensus_median_lfc_median < 0 & consensus_median_lfc_cutoff < 0 ~ TRUE,
      TRUE ~ FALSE
    )
  )

#save
write.csv(comparison_data_both, 
          paste0(output_dir, "/comparison_data_all_genes.csv"),
          row.names = FALSE)

#check differences in each direction
for (dir in c('up', 'down')) {
  
  comparison_data <- comparison_data_both %>%
    filter(final_direction == dir)
  
  # Genes in both for correlation plots
  genes_in_both <- comparison_data %>% filter(in_both)
  
  n_total <- nrow(comparison_data)
  n_both <- nrow(genes_in_both)
  n_only_median <- sum(comparison_data$only_median)
  n_only_cutoff <- sum(comparison_data$only_cutoff)
  n_consistent <- sum(genes_in_both$consistent_lfc, na.rm = TRUE)
  n_sig_both <- sum(genes_in_both$sig_both)
  
  #calculate spearman correlation of lfcs
  cor_lfc <- cor(genes_in_both$consensus_median_lfc_median,
                 genes_in_both$consensus_median_lfc_cutoff,
                 use = "complete.obs", method = "pearson")
  
  #plot correlation of shared genes across methods
  p1 <- ggplot(genes_in_both, aes(x = consensus_median_lfc_median,
                                  y = consensus_median_lfc_cutoff)) +
    geom_point(aes(color = -log10(min_pval)), size = 1.5, alpha = 0.9) +
    # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
    geom_smooth(method = "lm", se = TRUE, color = "darkred", alpha = 0.2, linewidth = 0.8) +
    scale_color_viridis_c(name = "-log10(min FDR)", option = "viridis") +
    labs(
      title = paste0("LFC comparison: ", dir, "-regulated genes"),
      subtitle = paste0("n = ", n_both, " genes | r = ", round(cor_lfc, 3),
                        " (Pearson)"),
      x = "Consensus median LFC (median method)",
      y = "Consensus median LFC (cutoff method)"
    ) +
    theme_bw(base_size = 12) +
    theme(aspect.ratio = 1)
  
  ggsave(paste0(output_dir, "/correlation_", dir, "_no_label.pdf"), plot = p1, width = 8, height = 8)
  write.csv(genes_in_both, paste0(output_dir, "/FigS4f_correlation_", dir, "_shared_in_both.csv"))
  
  # Save direction-specific data
  write.csv(comparison_data,
            paste0(output_dir, "/comparison_data_", dir, ".csv"),
            row.names = FALSE)
  
} #close up/dn loop

## plot best genes for both methods for direct comparison --> dotplot
genes_shared <- comparison_data_both %>%
  filter(in_both == TRUE)

#top 25 up and top 25 down based on highest average logfoldchange
top_genes <- genes_shared %>%
  arrange(lfc_mean) %>%
  {
    bind_rows(
      filter(., lfc_mean < 0) %>% slice_head(n = 25),
      filter(., lfc_mean > 0) %>% slice_tail(n = 25)
    )
  } %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(
    sig_median = p_aw_fdr_median < 0.05,
    sig_cutoff = p_aw_fdr_cutoff < 0.05
  )

#dumbbell plot - filled only for significant genes
top_genes <- top_genes %>%
  mutate(
    sig_median = p_aw_fdr_median < 0.05,
    sig_cutoff = p_aw_fdr_cutoff < 0.05
  )

p2 <- ggplot(top_genes, aes(y = fct_reorder(gene, lfc_mean))) +
  geom_segment(
    aes(x = consensus_median_lfc_median,
        xend = consensus_median_lfc_cutoff,
        yend = fct_reorder(gene, lfc_mean),
        color = final_direction),
    linewidth = 0.6,
    alpha = 0.6
  ) +
  # Median method points
  geom_point(
    aes(x = consensus_median_lfc_median,
        size = -log10(p_aw_fdr_median),
        fill = ifelse(sig_median, "#e09f3e", "white")),
    shape = 21, color = "#e09f3e", stroke = 1,
    alpha = 1
  ) +
  # Cutoff method points
  geom_point(
    aes(x = consensus_median_lfc_cutoff,
        size = -log10(p_aw_fdr_cutoff),
        fill = ifelse(sig_cutoff, "#335c67", "white")),
    shape = 21, color = "#335c67", stroke = 1,
    alpha = 1
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  scale_fill_identity() +
  scale_color_manual(
    values = c("up" = "#D55E00", "down" = "#0072B2"),
    name = "Direction"
  ) +
  scale_size_continuous(
    range = c(2, 6),
    name = "-log10(FDR)"
  ) +
  labs(
    title = "Top up/down-regulated genes: Method comparison",
    subtitle = paste0("Top 25 up + Top 25 down | Filled = FDR < 0.05 | Hollow = n.s.\n",
                      "Yellow = Median | Dark = Cutoff"),
    x = "Consensus median LFC",
    y = "Gene"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(paste0(output_dir, "/FigS4g_top_genes_dumbbell.pdf"), plot = p2, width = 10, height = 14)
