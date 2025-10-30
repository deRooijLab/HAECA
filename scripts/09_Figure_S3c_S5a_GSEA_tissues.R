#Geneset enrichment analysis
library(dplyr)
library(tidyverse)
library(clusterProfiler)
library(data.table)
library(org.Hs.eg.db)
library(ggplot2)
library(ComplexHeatmap)
library(colorRamp2)
library(RobustRankAggreg)
library(qvalue)
library(tidyplots)
library(cowplot)
library(ggrepel)

#### Load data ####
work_dir <- file.path("~/HAECA/global_analysis/deseq2")
setwd(work_dir)
data_dir <- paste0("~/HAECA/data")
output_dir <- paste0("~/HAECA/global_analysis/GSEA")
if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}


#source function
source("~/HAECA/scripts/GSEA_Biomex_function.R")

#get functions
load_DGEA <- function(tissue) {
  file_name <- file.path(work_dir, paste0("DGEA_", tissue, ".csv"))
  if (!file.exists(file_name)) return(NULL)
  
  DGEA <- read.csv(file_name, sep = ",", header = TRUE, dec = ".")[, 2:8]
  names(DGEA)[names(DGEA) == "gene"] <- "entity"
  names(DGEA)[names(DGEA) == "log2FoldChange"] <- "n"
  DGEA <- DGEA[, c("entity", "n")]
  
  features <- bitr(DGEA$entity, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  idx <- which(DGEA$entity %in% features$SYMBOL)
  DGEA <- DGEA[idx,]
  
  return(DGEA)
}

perform_gsea <- function(DGEA, genesets, geneset_name) {
  
  gsea <- biomex.performCompetitiveSetEnrichment(DGEA, featureSets = genesets, minimumSetSize = 5,
                                                 randomSeed = 1234)
  gsea <- gsea %>%
    mutate(Gene_symbols = sapply(Enriched_genes, function(x) suppressMessages(convert_entrez_to_symbols(x))))
  
  return(gsea)
  
}

# Prepare individual gene sets from files
geneset_names <- list('metabolic','senescence','mechanoregulation','hallmark')

#load all genesets
sets <- lapply(geneset_names, function(name) {
  file <- paste0(data_dir, "_", name, "_genesets_converted.csv")
  geneset <- read.csv(file)
  split(geneset$Gene, geneset$GeneSet)
})

# Assign category names as list names
names(sets) <- geneset_names

# Load data
files <- list.files(work_dir, pattern = ".csv", full.names = FALSE)
tissues <- unique(str_extract(files, "(?<=DGEA_)[a-zA-Z]+"))

#prepare list for combining results
gsea_results_list <- list()

for (tissue in tissues) {

    DGEA <- load_DGEA(tissue)
    if (is.null(DGEA)) next

    for (geneset_name in names(sets)) {
      
      genesets <- sets[[geneset_name]]
      gsea_results <- perform_gsea(DGEA, genesets, geneset_name)
      
      if (!is.null(gsea_results$errorMessage) && gsea_results$errorMessage == "Not enough sets") {
        break
      }
      
      #add to combined results list
      gsea_results_list[[tissue]][[geneset_name]] <- gsea_results %>%
      mutate(Tissue = tissue, Category = geneset_name)
      
      #save results for each tissue
      output_path <- file.path(output_dir, paste0("GSEA_", geneset_name))
      write.csv(gsea_results, file = paste0(output_path, ".csv"), row.names = FALSE)

    }

}

####Figure S5a - loop through through each geneset category to create separate heatmaps####
#function for heatmap
generate_heatmap <- function(data, geneset) {
  # Add significance column
  data <- data %>%
    mutate(Signif = case_when(
      Adjusted.pvalue < 0.001 ~ "***",
      Adjusted.pvalue < 0.01  ~ "**",
      Adjusted.pvalue < 0.05  ~ "*",
      TRUE                    ~ ""
    ))
  
  # Prepare heatmap data (NES values)
  heatmap_data <- data %>%
    group_by(Set, Tissue) %>%
    summarize(NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Tissue, values_from = NES, values_fill = NA)
  
  heatmap_matrix <- heatmap_data %>%
    arrange(Set) %>%
    column_to_rownames(var = "Set") %>%
    as.matrix() %>%
    t()
  
  # Prepare annotation data (significance asterisks)
  annot_data <- data %>%
    dplyr::select(Set, Tissue, Signif) %>%
    pivot_wider(names_from = Tissue, values_from = Signif, values_fill = "")
  
  #convert annotation data to a matrix with Set as rownames
  annot_text <- annot_data %>%
    arrange(Set) %>%
    column_to_rownames(var = "Set") %>%
    as.matrix() %>%
    t()
  
  #create heatmap
  heatmap <- Heatmap(
    heatmap_matrix,
    name = "NES",
    col = colorRamp2(c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
                     c("#296187", "white", "darkorange")),
    width = ncol(heatmap_matrix)*unit(5, "mm"),
    height = nrow(heatmap_matrix)*unit(5, "mm"),
    row_names_side = "left",
    column_names_side = "bottom",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(annot_text[i, j], x, y, gp = gpar(fontsize = 10, col = "black"))
    },
    column_title = paste(geneset, "genesets across tissues")
  )
  
  # Save to PDF
  pdf(file = paste0(output_dir, "/heatmap_", geneset, ".pdf"), width = 15, height = 8)
  draw(heatmap)
  dev.off()
}

#reformat for heatmap
combined_list <- list()

for (tissue in names(gsea_results_list)) {
  tissue_data <- gsea_results_list[[tissue]]
  tissue_combined <- bind_rows(
    lapply(names(tissue_data), function(category) {
      tissue_data[[category]] %>%
        mutate(Category = category)
    })
  ) %>%
    dplyr::mutate(Tissue = tissue)
  combined_list[[tissue]] <- tissue_combined
}

combined_genelist <- bind_rows(combined_list)

#Figure S5a - make heatmap for each geneset category
for (geneset in geneset_names) {
  category_data <- combined_genelist %>%
    filter(Category == geneset)

  #select sets significant in at least one tissue
  significant_sets <- category_data %>%
    group_by(Set) %>%
    filter(any(Adjusted.pvalue < 0.05)) %>%
    pull(Set) %>%
    unique()
  
  filtered_sig <- category_data %>%
    filter(Set %in% significant_sets)
  
  generate_heatmap(filtered_sig, geneset) #Figure S5a is Senescence heatmap
  
}

#Robust rank aggregation for all genesets by tissue
for(i in c("Up", "Down")) {
  if(i == "Up"){
ranked_lists <- combined_genelist %>%
  filter(NES > 0) %>%
  group_by(Tissue) %>%
  arrange(desc(NES)) %>% #sort from highest NES to lowest for increase in age
  mutate(Rank = row_number()) %>%
  dplyr::select(Tissue, Set, NES, Rank)
  } else {
    ranked_lists <- combined_genelist %>%
      filter(NES < 0) %>%
      group_by(Tissue) %>%
      arrange(NES) %>% #sort from lowest NES to highest for decrease in age
      mutate(Rank = row_number()) %>%
      dplyr::select(Tissue, Set, NES, Rank)
  }
  
# Create a named list where each tissue contains a rank
ranked_gene_sets <- ranked_lists %>%
  group_by(Tissue) %>%
  summarise(RankedGeneSets = list(Set)) %>%
  deframe()

#perform Robust Rank Aggregation
rra_results <- aggregateRanks(glist = ranked_gene_sets, method = "RRA")

# Adjust p-values using Benjamini-Hochberg (BH) method
rra_results$AdjustedPval <- p.adjust(rra_results$Score, method = "BH")

# Save results
write.csv(rra_results, file = paste0(output_dir, "/RRA_results_GSEA_tissues_", i, "regulated.csv"), row.names = FALSE)

}

#Figure S3c - make ranked plot for selected genesets after RRA
#set tissue colors
tissue_col <- c(
  'uterus' = "#d18c9f",
  'brain' = "#4FAE62",
  'lung' = "#F6C54D",
  'breast' = "#E37D46",
  'heart' = "#C02D45",
  'muscle' = "#8ecae6",
  'ovarium' = "#219ebc",
  'liver' = "#023047",
  'kidney' = "#dca100",
  'skin' = "#482374FF",
  'adipose' = "#fb8500",
  'intestine' = "#24878EFF",
  'bladder' = "#FF3F3F")

for(i in c("Up", "Down")) {
#get genesets that are significantly up/dn across n tissues
  shared_sets_filtered <- combined_genelist %>%
    left_join(y = shared_sets, by = "Set") %>%
    filter(
      n != '',
      Adjusted.pvalue < 0.05,
      if (i == "Up") NES > 0 else NES < 0
    ) %>%
    dplyr::select(Set, NES, Direction, Adjusted.pvalue, Tissue, n) %>%
    mutate(Labels = paste0(Set, " (", Tissue, ")")) 
    
#find top 3 genesets per tissue group
  if (i == "Down") {
    top_labels <- shared_sets_filtered %>%
    arrange(NES) %>%
      group_by(n) %>%
      do(head(., n = 3))
    
  } else {
    top_labels <- shared_sets_filtered %>%
      arrange(desc(NES)) %>%
      group_by(n) %>%
      do(head(., n = 3))
  }
  
  plot_shared <- ggplot(shared_sets_filtered, aes(x = factor(n), y = NES, color = factor(n))) +
    geom_point(size = 5, alpha = 0.8) +
    geom_label_repel(data = top_labels, aes(label = Labels),
                     color = "black", fill= NA, box.padding = 1.5) +
    {if(i == "Down") scale_y_reverse()} +
    scale_color_manual(values = rev(color_gradient)) +
    labs(
      title = paste0(i, "-regulated Significant Genesets Shared Across Tissues"),
      x = "Nr of Shared Tissues",
      y = "Normalized Enrichment Score"
    ) +
    theme_minimal() +
    theme(legend.position="none",
          axis.line = element_line(linewidth = 0.6, colour = "black", linetype=1),
          axis.text = element_text(color="black", size=12, face=1)
          ) +
    geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5), linetype = 2) +
    coord_flip()
  
  #save data from plot
  write.csv(shared_sets_filtered, file = paste0(output_dir, "/common_genesets_significant_", i, "regulated.csv"),
            row.names = FALSE, quote = FALSE)

  #save shared geneset plot
  pdf(file = paste0(output_dir, "/Figure_S3c_common_genesets_significant_", i, ".pdf"), width = 10, height = 10)
  print(plot_shared)
  dev.off()
  
}