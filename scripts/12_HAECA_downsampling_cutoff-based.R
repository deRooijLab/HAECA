library(tidyverse)
library(Seurat)
library(DESeq2)
library(edgeR)
library(future)
library(future.apply)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(ggplot2)
library(readr)
library(tidyr)
library(ggrepel)

#### Set directories
work_dir <- file.path("~/HAECA/global_analysis/tissues/filtered") #directory with all tissues as .RDS
setwd(work_dir)
output_dir <- file.path("~HAECA/global_analysis/downsampling/output_cutoff_parallel")
if(!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}

### MAIN SETTINGS
n_repeats <- 1000
output_suffix <- "age_brackets"

files <- list.files(work_dir, pattern = "_ECs_filtered.RDS", full.names = TRUE)
tissues <- basename(files) %>% str_remove("_ECs_filtered.RDS")

### 1: Build DESeq2 object from Seurat bulk
dds_from_bulk <- function(bulk_seurat, design, tissue) {
  
  metadata <- bulk_seurat@meta.data
  counts <- GetAssayData(bulk_seurat, slot = "counts", assay = "RNA")
  
  # check age groups
  tab <- table(metadata$Age_brackets)
  if (!all(c("young","aged") %in% names(tab))) {
    message("Skipping ", tissue, ": missing one age group")
    return(NULL)
  }
  if (any(tab[c("young","aged")] < 3)) {
    message("Skipping ", tissue, ": <3 samples per age group")
    return(NULL)
  }
  
  # filter low genes (CPM >10 in at least 10% samples)
  cpm_mat <- edgeR::cpm(counts)
  keep <- rowSums(cpm_mat > 10) >= 0.10 * ncol(cpm_mat)
  counts_filtered <- counts[keep, , drop=FALSE]
  
  # ensure no NA in design columns
    if (tissue == "kidney") {
    metadata <- metadata %>% drop_na(Age_brackets)
  } else {
    metadata <- metadata %>% drop_na(Age_brackets, Sampling_Study)
  }
  
  metadata$Age_brackets <- factor(metadata$Age_brackets, levels = c("young", "aged"))
  
  # create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = counts_filtered,
    colData   = metadata,
    design    = design
  )
  dds
}

### 2: subsample equal #cells per donor
sample_cells_per_donor <- function(seurat_obj, n_cells_per_donor, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  donors <- unique(seurat_obj$Donor)
  
  selected <- unlist(lapply(donors, function(d){
    cells <- colnames(seurat_obj)[seurat_obj$Donor == d]
    sample(cells, min(length(cells), n_cells_per_donor))
  }))
  subset(seurat_obj, cells = selected)
}

### PARALLEL PROCESSING SETUP
plan(multisession, workers = 12)

future_lapply(tissues, function(tissue) {
  
  tryCatch({
    
  message("\n### Processing tissue: ", tissue)
  
  # Load Seurat object
  seurat <- readRDS(file.path(work_dir, paste0(tissue, "_ECs_filtered.RDS")))
  seurat$Sampling_Study <- paste0(seurat$Sampling, "_", seurat$HAECA_ID)
  
  ### Define grouping
  grouping_vars <- if (tissue == "kidney") {
    c("Donor", "Age_brackets")
  } else {
    c("Donor", "Age_brackets", "Sampling_Study")
  }
  
  ### Define design
  design <- if (tissue == "kidney") {
    ~ Age_brackets
  } else {
    ~ Sampling_Study + Age_brackets
  }
  
  # Define age brackets: young = anything of 50 y/o or under, aged is 60 y/o or older
  seurat$Age_numeric <- as.numeric(seurat$Age_numeric)
  seurat@meta.data <- seurat@meta.data %>%
    mutate(Age_brackets = dplyr::case_when(
        seurat$Age_numeric <= 50 ~ "young",
        seurat$Age_numeric >= 60 ~ "aged",
        TRUE ~ "middle"
      )
    )

  # Full table of donors + numeric ages + bracket assignment
  Idents(seurat) <- "Age_brackets"
  table <- table(Idents(seurat), seurat$Age_numeric)
  write.csv(table, file = paste0(output_dir, "/Age_Donor_table_", tissue, "_subsampling.csv"))
  
  table <- table(Idents(seurat), seurat$HAECA_ID)
  write.csv(table, file = paste0(output_dir, "/Age_HAECA_ID_table_", tissue, "_subsampling.csv"))
  
  # drop the middle age bracket
  seurat <- subset(seurat, subset = Age_brackets %in% c("young","aged"))

  # Filter groups with enough cells
  donor_cell_count <- seurat@meta.data %>%
    tibble::as_tibble() %>%
    dplyr::select(all_of(grouping_vars)) %>%
    mutate(Group = apply(., 1, paste, collapse = "_")) %>%
    dplyr::count(Donor, Group, name = "CellCount")
  
  filtered_groups <- donor_cell_count %>%
    filter(CellCount >= 10) %>%
    pull(Group)

  seurat$Group <- apply(seurat@meta.data[, grouping_vars], 1, paste, collapse="_")
  seurat <- subset(seurat, Group %in% filtered_groups)

  # check age
  if (length(unique(seurat$Age_brackets)) < 2) {
    message("Skipping ", tissue, ": only one age group remains")
    return(NULL)
  }
  
  ### Determine subsampling depth (stable)
  target_cells_per_donor <- max(
    10,
    floor(quantile(donor_cell_count$CellCount, 0.25))
  )
  
  ### 1. ORIGINAL bulk analysis without subsampling
  bulk_real <- AggregateExpression(
    seurat, assays = "RNA", return.seurat = TRUE, group.by = grouping_vars
  )
  
  dds_real <- dds_from_bulk(bulk_real, design, tissue)
  if (is.null(dds_real)) return(NULL)
  
  dds_real <- DESeq(dds_real)
  real_results <- results(dds_real, contrast = c("Age_brackets","aged","young")) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    filter(!is.na(pvalue))
  
  write.csv(real_results,
            file.path(output_dir, paste0("DESeq2_", tissue, "_real.csv")),
            row.names = FALSE)
  
  ### 2. Repeated subsampling
  message("Running ", n_repeats, " subsampling runs...")
  
  repeat_results <- vector("list", n_repeats)
  
  for (rep in seq_len(n_repeats)) {
    seurat_rep <- sample_cells_per_donor(seurat, target_cells_per_donor, seed = rep)
    
    bulk_rep <- AggregateExpression(
      seurat_rep, assays="RNA", return.seurat=TRUE, group.by=grouping_vars
    )
    
    dds_rep <- dds_from_bulk(bulk_rep, design, tissue)
    if (is.null(dds_rep)) next
    
    dds_rep <- DESeq(dds_rep)
    res_rep <- results(dds_rep, contrast=c("Age_brackets","aged","young")) %>%
      as.data.frame() %>% rownames_to_column("gene") %>% filter(!is.na(pvalue))
    
    repeat_results[[rep]] <- res_rep
  }
  
  repeat_all <- bind_rows(lapply(seq_along(repeat_results), function(i){
    if (is.null(repeat_results[[i]])) return(NULL)
    repeat_results[[i]] %>%
      mutate(
        Tissue = tissue,
        repeat_id = i,
        run_type = "repeat"
      )
  }))
  
  write.csv(repeat_all,
            file.path(output_dir, paste0("DESeq2_", tissue, "_allRepeats.csv")),
            row.names=FALSE)
  
  ### 3. Stability metrics
  gene_stability <- repeat_all %>%
    group_by(Tissue, gene) %>%
    summarise(
      times_significant = sum(pvalue < 0.25, na.rm=TRUE),
      prop_sig = times_significant / n_repeats,
      mean_lfc = mean(log2FoldChange, na.rm=TRUE),
      sd_lfc = sd(log2FoldChange, na.rm=TRUE),
      median_p = median(pvalue, na.rm=TRUE)
    ) %>% ungroup()
  
  write.csv(gene_stability,
            file.path(output_dir, paste0("DESeq2_gene_stability_", tissue, "_subsampling.csv")),
            row.names=FALSE)

  # Save all repeat results
  write.csv(repeat_all, paste0(output_dir, "/DESeq2_", tissue, "_allRepeats.csv"), row.names = FALSE)
  
  # GeneStats: Summarize repeated subsampling results per gene
  gene_stats <- repeat_all %>%
    group_by(gene) %>%
    summarise(
      prop_sig = mean(pvalue < 0.25, na.rm = TRUE),
      mean_lfc = mean(log2FoldChange, na.rm = TRUE),
      sd_lfc   = sd(log2FoldChange, na.rm = TRUE),
      median_p = median(pvalue, na.rm = TRUE),
      median_padj = median(padj, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    arrange(median_p) %>%
    mutate(
      gene_order = row_number(),
      neglog10_p = -log10(median_p)
    )
  
  write.csv(gene_stats, paste0(output_dir, "/GeneStats_", tissue, "_allRepeats.csv"), row.names = FALSE)
  
  message("Finished ", tissue)
  
  return(TRUE)
  
  }, error = function(e) {
    
    message("Error in tissue ", tissue , ": ", conditionMessage(e))
    message("Skipping and continuing with next tissue.")
    
    return(NULL)   # or return(list(error=t))
  })
  
})
