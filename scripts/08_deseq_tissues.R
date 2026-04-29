#Wald test (DESeq2) for each tissue and age bracket comparison

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