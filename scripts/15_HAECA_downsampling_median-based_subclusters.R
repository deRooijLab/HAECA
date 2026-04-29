library(tidyverse)
library(Seurat)
library(DESeq2)
library(edgeR)
library(future)
library(future.apply)
library(ggrepel)
library(ggplot2)

#### Set directories
work_dir <- file.path("~/HAECA/global_analysis/tissues/filtered")
setwd(work_dir)
output_dir <- file.path("~HAECA/global_analysis/downsampling_tissue_subclusters/output_median")
dir.create(output_dir,recursive = T, showWarnings = T)

### MAIN SETTINGS
n_repeats <- 1000
output_suffix <- "age_brackets"

work_files <- list.files(work_dir, pattern = "_ECs_filtered.RDS", full.names = TRUE)
tissues <- basename(work_files) %>% str_remove("_ECs_filtered.RDS")

### 1: Build DESeq2 object from Seurat bulk
dds_from_bulk <- function(bulk_seurat, design, tissue, cluster) {
  
  metadata <- bulk_seurat@meta.data
  counts <- GetAssayData(bulk_seurat, layer = "counts", assay = "RNA")
  
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
  
  dds$Age_brackets <- relevel(dds$Age_brackets, ref = "young")
  
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

### 3: Count cells and donors for reporting
count_cells_and_donors <- function(seurat_obj, tissue, cluster, grouping_vars) {
  
  meta <- seurat_obj@meta.data
  
  # Donor counts per age bracket
  donor_summary <- meta %>%
    group_by(Age_brackets) %>%
    summarise(
      n_donors = n_distinct(Donor),
      n_cells_total = n(),
      .groups = "drop"
    ) %>%
    mutate(
      tissue = tissue,
      cluster = cluster
    )
  
  # Cell counts per donor per group
  cell_summary <- meta %>%
    mutate(Group = apply(select(., all_of(grouping_vars)), 1, paste, collapse = "_")) %>%
    group_by(Donor, Group, Age_brackets) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    mutate(
      tissue = tissue,
      cluster = cluster
    )
  
  # Detailed summary for QC
  qc_summary <- meta %>%
    group_by(Age_brackets, Donor) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(Age_brackets) %>%
    summarise(
      n_donors = n(),
      total_cells = sum(n_cells),
      mean_cells_per_donor = mean(n_cells),
      median_cells_per_donor = median(n_cells),
      min_cells_per_donor = min(n_cells),
      max_cells_per_donor = max(n_cells),
      .groups = "drop"
    ) %>%
    mutate(
      tissue = tissue,
      cluster = cluster
    )
  
  list(
    donor_summary = donor_summary,
    cell_summary = cell_summary,
    qc_summary = qc_summary
  )
}

### prepare dataframes
all_donor_summaries <- data.frame()
all_cell_summaries <- data.frame()
all_qc_summaries <- data.frame()

all_run_logs <- tibble::tibble(
  tissue = character(),
  cluster = character(),
  n_cells_total = integer(),
  n_donors_total = integer(),
  n_donors_young = integer(),
  n_donors_aged = integer(),
  n_cells_young = integer(),
  n_cells_aged = integer(),
  target_cells_per_donor = integer(),
  n_successful_repeats = integer(),
  n_genes_tested = integer(),
  n_genes_sig_real = integer(),
  status = character()
)

### PARALLEL PROCESSING SETUP
plan(multisession, workers = 12)

results_list <- future_lapply(tissues, function(tissue) {
  
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
  
  # Get available clusters for this tissue
  available_clusters <- levels(seurat$clustering_global)
  available_clusters <- available_clusters[available_clusters %in% unique(seurat$clustering_global)]
  
  # Results containers for this tissue
  tissue_results <- list()
  
  #subset for each cluster
  for(cluster in available_clusters) {
    
    tryCatch({
      
    print(paste0("now starting bulk ", tissue, " ", cluster))
    
    i <- subset(x = seurat, subset = clustering_global == cluster)
    
    # Skip if too few cells
    if (ncol(i) < 30) {
      message("    Skipping ", tissue, "_", cluster, ": only ", ncol(i), " cells")
      next
    }
    
    counts_info <- count_cells_and_donors(i, tissue, cluster, grouping_vars)
    
    # Check minimum requirements
    donor_check <- counts_info$donor_summary
    if (nrow(donor_check) < 2 || any(donor_check$n_donors < 3)) {
      message("    Skipping ", tissue, "_", cluster, ": insufficient donors per age group")
      message("    Donors: ", paste(donor_check$Age_brackets, "=", donor_check$n_donors, collapse = ", "))
      next
    }
  
    donor_cell_df <- i@meta.data %>%
      mutate(Group = apply(select(., all_of(grouping_vars)), 1, paste, collapse = "_")) %>%
      dplyr::count(Donor, Group, name = "CellCount")
    
    valid_donors <- donor_cell_df %>% 
      filter(CellCount >= 10) %>% 
      pull(Donor) %>% 
      unique()
    
    # Subset to valid donors
    i <- subset(i, Donor %in% valid_donors)
    
    # Re-check after filtering
    if (ncol(i) < 50) {
      message("    Skipping ", tissue, "_", cluster, ": only ", ncol(i), " cells after donor filter")
      next
    }
    
    # Re-count after filtering
    counts_info_filtered <- count_cells_and_donors(i, tissue, cluster, grouping_vars)
    donor_check_filtered <- counts_info_filtered$donor_summary
    
    if (nrow(donor_check_filtered) < 2 || any(donor_check_filtered$n_donors < 3)) {
      message("    Skipping ", tissue, "_", cluster, ": insufficient donors after filtering")
      next
    }
    
    message("    ✓ Passed QC: ", sum(donor_check_filtered$n_cells_total), " cells, ",
            sum(donor_check_filtered$n_donors), " donors")
  
  ### Determine subsampling depth (stable)
    cells_per_donor <- donor_cell_df %>%
      filter(Donor %in% valid_donors) %>%
      group_by(Donor) %>%
      summarise(total_cells = sum(CellCount), .groups = "drop")
    
  target_cells_per_donor <- max(
    10,
    floor(quantile(donor_cell_df$CellCount, 0.25))
  )
  
  message("    Target cells per donor for subsampling: ", target_cells_per_donor)
  
  ### 1. ORIGINAL bulk analysis without subsampling
  bulk_real <- AggregateExpression(
    i, assays = "RNA", return.seurat = TRUE, group.by = grouping_vars
  )
  
  dds_real <- dds_from_bulk(bulk_real, design, tissue, cluster)
  if (is.null(dds_real)) {
    message("    Skipping ", tissue, "_", cluster, ": DESeq2 object creation failed")
    next
  }
  
  dds_real <- DESeq(dds_real)
  real_results <- results(dds_real, contrast = c("Age_brackets", "aged", "young")) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    filter(!is.na(pvalue)) %>%
    mutate(tissue = tissue,
      cluster = cluster)
  
  # Save real results
  write.csv(real_results,
            file.path(output_dir, paste0("DESeq2_", tissue, "_", cluster, "_real.csv")),
            row.names = FALSE, quote = FALSE)
  
  ### 2. Repeated subsampling
  message("Running ", n_repeats, " subsampling runs...")
  
  repeat_results <- vector("list", n_repeats)
  
  for (rep in seq_len(n_repeats)) {
    
    seurat_rep <- sample_cells_per_donor(i, target_cells_per_donor, seed = rep)
    
    # Skip if too few cells after subsampling
    if (ncol(seurat_rep) < 15) next
    
    bulk_rep <- AggregateExpression(
      seurat_rep, assays="RNA", return.seurat=TRUE, group.by=grouping_vars
    )
    
    dds_rep <- dds_from_bulk(bulk_rep, design, tissue, cluster)
    if (is.null(dds_rep)) next
    
    dds_rep <- DESeq(dds_rep, quiet = TRUE)
    res_rep <- results(dds_rep, contrast=c("Age_brackets","aged","young")) %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      filter(!is.na(pvalue))
    
    repeat_results[[rep]] <- res_rep
    
    if (rep %% 100 == 0) message("      Completed ", rep, "/", n_repeats)
  }
  
  repeat_all <- bind_rows(lapply(seq_along(repeat_results), function(r){
    if (is.null(repeat_results[[r]])) return(NULL)
    repeat_results[[r]] %>%
      mutate(
        Tissue = tissue,
        cluster = cluster,
        repeat_id = r,
        run_type = "repeat"
      )
  }))
  
  if (nrow(repeat_all) == 0) {
    message("    Warning: No successful subsampling runs for ", tissue, "_", cluster)
    next
  }
  
  write.csv(repeat_all,
            file.path(output_dir, paste0("DESeq2_", tissue, "_", cluster, "_allRepeats.csv")),
            row.names=FALSE)
  
  ### 3. Stability metrics
  gene_stability <- repeat_all %>%
    group_by(gene) %>%
    summarise(
      times_significant = sum(pvalue < 0.25, na.rm=TRUE),
      prop_sig = times_significant / n_repeats,
      mean_lfc = mean(log2FoldChange, na.rm=TRUE),
      sd_lfc = sd(log2FoldChange, na.rm=TRUE),
      median_p = median(pvalue, na.rm=TRUE),
      median_padj = median(padj, na.rm = TRUE)
    ) %>%
    mutate(
      tissue = tissue,
      cluster = cluster) %>%
    ungroup()
  
  write.csv(gene_stability,
            file.path(output_dir, paste0("DESeq2_gene_stability_", tissue, "_", cluster, "_subsampling.csv")),
            row.names=FALSE)

  ### 4. final report
  tissue_results[[cluster]] <- list(
    donor_summary = counts_info_filtered$donor_summary,
    cell_summary = counts_info_filtered$cell_summary,
    qc_summary = counts_info_filtered$qc_summary,
    run_log = data.frame(
      tissue = tissue,
      cluster = cluster,
      n_cells_total = ncol(i),
      n_donors_total = length(valid_donors),
      n_donors_young = donor_check_filtered$n_donors[donor_check_filtered$Age_brackets == "young"],
      n_donors_aged = donor_check_filtered$n_donors[donor_check_filtered$Age_brackets == "aged"],
      n_cells_young = donor_check_filtered$n_cells_total[donor_check_filtered$Age_brackets == "young"],
      n_cells_aged = donor_check_filtered$n_cells_total[donor_check_filtered$Age_brackets == "aged"],
      target_cells_per_donor = target_cells_per_donor,
      n_successful_repeats = sum(!sapply(repeat_results, is.null)),
      n_genes_tested = nrow(real_results),
      n_genes_sig_real = sum(real_results$pvalue < 0.05, na.rm = TRUE),
      status = "completed"
    )
  )
  
  message("Completed: ", tissue, "_", cluster)
  
    }, error = function(e) {
      message("Error in ", tissue, "_", cluster, ": ", conditionMessage(e))
      next
    })
  }  # End cluster loop
  
  return(tissue_results)
  
  }, error = function(e) {
    message("Error in tissue ", tissue, ": ", conditionMessage(e))
    return(NULL)
  })
  
}, future.seed = TRUE, future.stdout = TRUE, future.conditions = "message")