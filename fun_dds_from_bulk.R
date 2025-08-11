#function for dds from bulk
dds_from_bulk <- function(tissue, design) {
  
  file_name <- paste0("bulk_", tissue, "_age_brackets.rds")
  bulk <- readRDS(file_name)
  
  metadata <- bulk@meta.data
  
  counts <- GetAssayData(bulk, assay  ="RNA", layer = 'counts')
  
  # Only proceed if every age bracket has at least 4 donors (for stats)
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
  
  saveRDS(dds, paste0(output_dir, "/dds_",tissue ,".RDS"))
  
  return(dds)
  
}
