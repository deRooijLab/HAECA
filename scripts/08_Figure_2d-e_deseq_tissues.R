#Wald test for each tissue and age bracket comparison
#combined RRA analysis for Figures 2d and 2e

library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(tidyverse)
library(edgeR)
library(DESeq2)
library(EnhancedVolcano)
library(cowplot)
library(RColorBrewer)
library(RobustRankAggreg)
library(qvalue)
library(openxlsx)

#### Load data ####
work_dir <- file.path("~/HAECA/global_analysis/bulk")
setwd(work_dir)
output_dir <- paste0("~/HAECA/global_analysis/deseq2")
if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
rra_dir <- file.path(paste0(output_dir, "/RRA_combined"))
if(!dir.exists(rra_dir)){dir.create(rra_dir,recursive = T)}

#function for dds from bulk
dds_from_bulk <- function(tissue, design) {
  
  file_name <- paste0("bulk_", tissue, "rds")
  bulk <- readRDS(file_name)
  
  metadata <- bulk@meta.data
  
  counts <- GetAssayData(bulk, assay  ="RNA", layer = 'counts')
  
  # Only proceed if every age bracket has at least 4 donors, for stats
  age_counts <- table(bulk$Age_brackets)
  
  #report error if less than 4 donors in each age bracket
  if (!all(c("aged", "young") %in% names(age_counts)) ||
      any(age_counts[c("aged", "young")] < 4)) {
    message(paste("Insufficient samples for tissue:", tissue))
    return(NULL)
  }
  
  # Filter out lowly expressed genes
  myCPM <- cpm(counts)
  thresh <- myCPM > 10
  to_keep <- 0.10*ncol(myCPM) #keep genes detected above thresh in 10% of samples
  keep <- rowSums(thresh) >= to_keep
  
  # Subset the rows of countdata to keep the more highly expressed genes
  counts_filtered <- counts[keep,]
  dim(counts_filtered)
  
  #perform differential expression analysis, comparing the different age brackets
  dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = metadata,
                                design = design)
  
  dds$Age_brackets <- relevel(dds$Age_brackets, ref = "young")
  
  saveRDS(dds, paste0(output_dir, "/dds_", tissue, ".RDS"))
  
  return(dds)
  
}

cols <- c("#266c91","#f0cf16")

# Load data
files <- list.files(work_dir, pattern = ".rds", full.names = FALSE)
tissues <- unique(str_extract(files, "(?<=bulk_)[a-zA-Z]+"))

contrast <- c("Age_brackets", "aged", "young")

#prepare list for collection of tissue results
combined_list <- list()

#deseq for each tissue
for(tissue in tissues) {
  
  if (tissue == "kidney") {
    design <- as.formula(~ Age_brackets)
    reduced <- as.formula(~ 1)
  } else {
    design <- as.formula(~ Sampling_Study + Age_brackets)
    reduced <- as.formula(~ Sampling_Study)
  }
  
  dds <- dds_from_bulk(tissue = tissue, design = design)
  
  if (is.null(dds)) {
    message(paste("Skipping", tissue))
  } else {
  
    message(paste("Running Wald test for", tissue))
    dds_wald <- DESeq(dds)
    results <- results(dds_wald, contrast = contrast)

  results_tb <- results %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  #add to list
  combined_list[[tissue]] <- results_tb
  
  #save results
  write.csv(results_tb, file = paste0(output_dir, "/DGEA_", tissue, ".csv"))
  }
}

###Robust Rank Aggregation from DGEA####

#define colors for plot
color_gradient <- rev(c("#D1422F","#F4AB5C","#e9c46a","#ACC8BE","#1A5B5B","#264653"))

tissue_nr_all <- list()

for (direction in c("up", "down")) {
  
  filtered_genes <- filter_genes(data = combined_list, logfc_threshold = 0, direction = direction)
  
  # Prepare wide format
  gene_df_wide <- create_wide_gene_df(filtered_genes)
  
  # Extract gene sets
  gene_sets <- lapply(filtered_genes, function(x) x$gene)
  names(gene_sets) <- names(filtered_genes)
  
  # Create gene matrix
  gene_matrix <- create_gene_matrix(gene_sets)
  
  # Perform RRA analysis
  glist <- setNames(lapply(filtered_genes, function(x) x$gene), names(filtered_genes))
  ar <- perform_rra_analysis(glist)
  
  #pick top genes
  top_genes <- ar %>%
    top_n(n=25, -qvalue) %>%
    mutate(qvalue_log = -log10(qvalue))
  
  # Generate plots
  plot <- plot_point_chart(top_genes, color_gradient)
  
  tissue_nr <- gene_matrix %>%
    mutate(n = rowSums(.)) %>%
    rownames_to_column(var = 'Gene') %>%
    left_join(ar %>% dplyr::select(Name, qvalue) %>% mutate(qvalue = -log10(qvalue)), by = c("Gene" = "Name"))
  
  #create nr of tissues per gene for table
  tissue_nr_all[[direction]] <- tissue_nr
  
  #Table S6 - save RRA results
  write.csv(tissue_nr, file = paste0(rra_dir, '/RRA_', direction, '.csv'), row.names = F, quote = F)
  
  # Save plots
  pdf(file = paste0(rra_dir, '/Figure_2d-e_', direction, 'regulated.pdf'), width = 10, height = 10)
  print(plot)
  dev.off()
  
}

#save together as excel sheet
wb <- createWorkbook()
for (name in names(tissue_nr_all)) {
  df <- tissue_nr_all[[name]]
    #add to excelsheets
  addWorksheet(wb, name)
  writeData(wb, name, df)
}
saveWorkbook(wb, file.path(paste0(rra_dir, "/wald_RRA_meta_analysis.xlsx")), overwrite = T)

