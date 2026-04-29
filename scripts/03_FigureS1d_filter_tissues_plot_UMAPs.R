library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(data.table)
library(tidyverse)
library(viridis)
library(presto)
library(openxlsx)

#set directory where .RDS files should be changed
work_dir <- file.path("~/HAECA/global_analysis/tissues")
setwd(work_dir)
output_dir <- file.path(paste0(work_dir, "/filtered"))
if(!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
marker_dir <- file.path(paste0(output_dir, "/markers"))
if(!dir.exists(marker_dir)){dir.create(marker_dir, recursive = T)}
plot_dir <- file.path(paste0(output_dir, "/visualization"))
if(!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

#load biomaRt annotations
ensembl <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                            dataset = "hsapiens_gene_ensembl",
                            host = "http://www.ensembl.org")

annotations <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"), mart = ensembl)

ambient_genes <- read.csv(paste0(script_dir, "/ambient_genes.csv"))[,1]

#list files in directory
files <- list.files(work_dir, recursive = TRUE, pattern = "annotated_prefinal.RDS")

#call functions
#function to add and calculate data age info
HAECA_metadata_formatting <- function(data) {
  
  #1. calculate median in case of given range
  metadata <- data@meta.data %>%
    mutate(Age_numeric = sapply(Age, convert_age_range_to_median)) %>%
    rownames_to_column(var = "cell_id")
  
  #2. calculate median across whole tissue depending on Donor ages
  median <- metadata %>%
    distinct(Donor, Age_numeric) %>%
    mutate(
      median_age = median(Age_numeric, na.rm = TRUE),
      cutoff_age = ifelse(median_age > 60, 60, median_age)
    ) %>%
    ungroup() %>%
    select(-Age_numeric)
  
  #3. add age brackets depending on median cutoff
  metadata <- metadata %>%
    left_join(median, by = "Donor") %>%
    mutate(Age_brackets = ifelse(Age_numeric < cutoff_age, "young", "aged")) %>%
    ungroup()
  
  #4. calculate decade range
  metadata <- metadata %>%
    mutate(Decade = paste0(10 * (Age_numeric %/% 10), "-", 10 * (Age_numeric %/% 10) + 9))
  
  #5. add new columns to data
  metadata <- metadata %>%
    column_to_rownames(var = "cell_id") %>%
    select(Age_brackets,Age_numeric,Decade,median_age,cutoff_age)
  
  ids.metadata <- row.names(metadata)
  ids.seurat <- row.names(data@meta.data)
  
  if(identical(ids.metadata, ids.seurat)){
    data <- AddMetaData(object = data, metadata = metadata)
    
  } else {
    stop("Error: rownames of seurat object and formatted metadata not identical.")
  }
  
  #6. reduce metadata on necessary info only
  colnames(data@meta.data)
  metadata_include <- c("orig.ident","nCount_RNA","nFeature_RNA","Genome_build","Donor","Age",
                        "Sampling","Sample","Method","Sample_info","Gender","Tissue",
                        "GEO_accession","HAECA_ID","percent.mt","seurat_clusters",
                        "Cell_source","unintegrated_clusters","clustering_fine","clustering_global",
                        "Age_brackets","Age_numeric","Decade","clustering_merged","Cell_type_authors",
                        "median_age", "cutoff_age"
  )
  
  data@meta.data <- data@meta.data[,names(data@meta.data) %in% metadata_include]
  data@meta.data <- droplevels(data@meta.data)
  
  return(data)
  
}

#make function to filter gene lists
HAECA_gene_filtering <- function(data, annotations, ambient_genes) {
  
  #1. filter to include only protein-coding genes
  #look for overlap in genes
  annotations_seurat <- intersect(as.character(rownames(data)), as.character(annotations$hgnc_symbol))
  annotations_filtered <- annotations %>% 
    filter(hgnc_symbol %in% annotations_seurat)
  
  #select protein_coding as gene_biotype of interest
  prot <- which(annotations_filtered$gene_biotype == "protein_coding")
  final_genes <- annotations_filtered[prot,]
  
  data <- subset(data, features = final_genes$hgnc_symbol, invert = FALSE)
  
  #2. filter ambient genes from curated list across all studies
  #remaining unwanted genes (antisense, micro, lnc, non-coding RNAs, ENSG annotations)
  #filter out all mitochondrial and ribosomal genes
  
  seurat_genes <- as.character(rownames(data))
  
  remove_ambient <- intersect(seurat_genes, ambient_genes)
  
  remove_pattern <- seurat_genes %>%
    as_tibble() %>%
    filter(str_detect(value, "^(RPS|RPL|MT-)|-AS1|ENSG0000|LINC|MIR|RF000")) %>%
    pull()
  
  #combine both
  remove_genes <- unique(c(remove_ambient, remove_pattern))
  remaining_genes <- setdiff(seurat_genes, remove_genes)
  
  data <- subset(data, features = remaining_genes, invert = FALSE)
  
  #3. filter to intersect across all studies (HAECA_ID)
  # Only work with genes detected in every HAECA_ID (to avoid study-specific expressions appearing in the analysis)
  Idents(data) <- "HAECA_ID"
  IDs <- data@meta.data %>%
    distinct(HAECA_ID) %>%
    pull() %>%
    as.vector()
  
  filtered_HAECA_ID <- lapply(IDs, function(id) {
    seurat_subset <- subset(data, idents = id)
    idx <- which(rowSums(GetAssayData(seurat_subset, layer = "counts")) > 0)
    seurat_subset[idx, ]
  })
  
  intersection <- Reduce(intersect, lapply(filtered_HAECA_ID, rownames))
  data <- subset(data, features = intersection)
  
  return(data)
  
}

#convert_age_range_to_median function 
convert_age_range_to_median <- function(age_range) {
  # Remove any non-numeric characters like "_y"
  age_range_clean <- gsub("[^0-9-]", "", age_range)
  
  if (grepl("-", age_range_clean)) {
    # If given as range, split and calculate the median
    age_split <- as.numeric(strsplit(age_range_clean, "-")[[1]])
    
    #check if valid
    if (length(age_split) == 2 && !any(is.na(age_split))) {
      return(mean(age_split))
    } else {
      warning(paste("Invalid age range:", age_range))
      return(NA)  # Return NA if the range is invalid
    }
  } else {
    #return the numeric value
    single_age <- as.numeric(age_range_clean)
    
    if (!is.na(single_age)) {
      return(single_age)
    } else {
      warning(paste("Invalid single age:", age_range))
      return(NA)  # Return NA if the value is invalid
    }
  }
}

#call tissue color mappings
tissue_cluster_colors <- function(){
  
  # Store results in a list
  tissue_color_mappings <- list()
  
  tissue_color_mappings[["adipose"]] <- c(
    "artery_1" = "#87092b",
    "artery_2" = "#f70230",
    "capillary_arterial" = "#d18c9f",
    "capillary_1" = "#8cd1a9",
    "capillary_2" = "#72a174",
    "capillary_3" = "#327878",
    "capillary_4" = "#bbeded",
    "vein_1" = "#269feb",
    "vein_2" = "#093c80",
    "vein_3" = "#8565a6",
    "lymphatic" = "#cfca42"
  )
  
  tissue_color_mappings[["brain"]] <- c(
    "artery_1" = "#87092b",
    "artery_2" = "#f70230",
    "artery_3" = "#D08B9F",
    "capillary_arterial" = "#f5dce0",
    "capillary_angiogenic" = "#D15C2C",
    "capillary_1" = "#8cd1a9",
    "capillary_2" = "#72a174",
    "capillary_venous_1" = "#16C4C1",
    "capillary_venous_2" = "#6fc3ed",
    "capillary_venous_3" = "lightgrey",
    "vein_1" = "#093c80",
    "vein_2" = "#0e5ec7",
    "vein_3" = "#8565a6",
    "vein_4" = "#269feb",
    "choroid_plexus" = "#682860"
  )
  
  tissue_color_mappings[["breast"]] <- c(
    "artery_1" = "#87092b",
    "artery_2" = "#f70230",
    "capillary_arterial_1" = "#D08B9F",
    "capillary_arterial_2" = "#f5dce0",
    "capillary_arterial_3" = "#c75077",
    "capillary_angiogenic" = "#D15C2C",
    "capillary_1" = "#8cd1a9",
    "capillary_2" = "#72a174",
    "vein_1" = "#093c80",
    "vein_2" = "#0e5ec7",
    "vein_3" = "#269feb",
    "vein_4" = "lightgrey",
    "lymphatic_1" = "#cfca42",
    "lymphatic_2" = "#dbac3d"
  )
  
  tissue_color_mappings[["heart"]] <- c(
    "artery_1" = "#87092b",
    "artery_2" = "#f70230",
    "capillary_arterial_1" = "#d18c9f",
    "capillary_arterial_2" = "#f5dce0",
    "capillary_endocardial" = "#86608E",
    "capillary_1" = "#8cd1a9",
    "capillary_2" = "#72a174",
    "capillary_3" = "#327878",
    "capillary_4" = "#bbeded",
    "capillary_5" = "#0D6153",
    "capillary_6" = "#7BA05B",
    "capillary_venous" = "#16C4C1",
    "vein_1" = "#269feb",
    "vein_2" = "#093c80",
    "lymphatic_1" = "#cfca42",
    "lymphatic_2" = "#dbac3d"
  )
  
  tissue_color_mappings[["intestine"]] <- c(
    "artery_1" ="#87092b",
    "artery_2" = "#f70230",
    "capillary_1" = "#8cd1a9", 
    "capillary_2" = "#72a174", 
    "capillary_3" = "#327878",
    "capillary_4" = "#bbeded",
    "vein_1" = "#269feb",
    "vein_2" = "#093c80",
    "vein_3" = "#8565a6",
    "lymphatic_1" = "#cfca42",
    "lymphatic_2" = "#dbac3d",
    "lymphatic_3" = "#FFEF00",
    "lymphatic_4" = "#cfca42"
  )
  
  tissue_color_mappings[["kidney"]] <- c(
    "artery" = "#87092b",
    "capillary_arterial_DVR_1" = "#d18c9f",
    "capillary_arterial_DVR_2" = "#f5dce0",
    "capillary_arterial_afferent/efferent" = "#F4C2C2",
    "capillary_angiogenic" = "#D15C2C",
    "capillary_peritubular_1" = "#8cd1a9",
    "capillary_peritubular_2" = "#72a174",
    "capillary_glomerular" = "#86608E",
    "capillary_venous_AVR_1" = "#16C4C1",
    "capillary_venous_AVR_2" = "#6fc3ed",
    "lymphatic" = "#cfca42"
  )
  
  tissue_color_mappings[["liver"]] <- c(
    "artery_large" = "#87092b",
    "artery_periportal" = "#d18c9f",
    "capillary_peribiliary" = "#8cd1a9",
    "capillary_LSEC_1" = "#86608E",
    "capillary_LSEC_2" = "#682860",
    "capillary_LSEC_3" = "#B284BE",
    "vein_periportal_1" = "#6fc3ed",
    "vein_periportal_2" = "#269feb",
    "vein_central" = "#093c80",
    "lymphatic" = "#cfca42"
  )
  
  tissue_color_mappings[["lung"]] <- c(
    "artery" = "#87092b",
    "capillary_arterial" = "#d18c9f",
    "capillary_angiogenic" = "#D15C2C",
    "capillary_aerocyte" = "#86608E",
    "capillary_general_1" = "#8cd1a9",
    "capillary_general_2" = "#72a174",
    "capillary_systemic" = "#682860",
    "vein_systemic" = "#6b8ab2",
    "vein_1" = "#269feb",
    "vein_2" = "#093c80",
    "lymphatic_1" = "#cfca42",
    "lymphatic_2" = "#dbac3d"
  )
  
  tissue_color_mappings[["muscle"]] <- c(
    "artery" = "#87092b",
    "capillary_arterial" = "#d18c9f",
    "capillary_angiogenic" = "#D15C2C",
    "capillary_1" = "#8cd1a9",
    "capillary_2" = "#72a174", 
    "vein_1" = "#269feb",
    "vein_2" = "#093c80",
    "vein_3" = "#8565a6",
    "vein_4" = "#6b8ab2",
    "vein_5" ="#0e5ec7",
    "lymphatic_1" = "#cfca42",
    "lymphatic_2" = "#dbac3d")
  
  tissue_color_mappings[["ovarium"]] <- c(
    "artery_1" ="#87092b",
    "artery_2" = "#f70230",
    "capillary_arterial" = "#d18c9f",
    "capillary_angiogenic" = "#D15C2C",
    "capillary" = "#8cd1a9",
    "capillary_venous_1" = "#16C4C1", 
    "capillary_venous_2" ="#6fc3ed", 
    "capillary_venous_3" ="#50c7c7",
    "vein" = "#269feb",
    "lymphatic" ="#cfca42")
  
  tissue_color_mappings[["skin"]] <- c(
    "artery" = "#87092b",
    "capillary_arterial" = "#d18c9f",
    "capillary_angiogenic" = "#D15C2C",
    "capillary_1" = "#8cd1a9",
    "capillary_2" = "#72a174",
    "capillary_venous" = "#16C4C1",
    "vein_1" = "#269feb",
    "vein_2" = "#093c80",
    "vein_3" = "#8565a6",
    "vein_4" = "#6b8ab2",
    "lymphatic_1" = "#cfca42",
    "lymphatic_2" = "#dbac3d"
  )
  
  tissue_color_mappings[["uterus"]] <- c(
    "artery" = "#87092b",
    "capillary_arterial_1" = "#d18c9f",
    "capillary_arterial_2" = "#f5dce0",
    "capillary_1" = "#8cd1a9",
    "capillary_2" = "#72a174",
    "capillary_venous_1" = "#16C4C1",
    "capillary_venous_2" = "#6fc3ed",
    "capillary_venous_3" = "#50c7c7",
    "vein_1" = "#269feb",
    "vein_2" = "#093c80",
    "vein_3" = "#8565a6",
    "vein_4" = "#0e5ec7"
  )
  
  return(tissue_color_mappings)
}

#call function to calculate markers
HAECA_top5_markers <- function(seurat) {
  
  #store results in list
  results <- list()
  
  #switch to RNA assay for marker calculation
  DefaultAssay(seurat) <- "RNA"
  
  #if clustering_merged present, here calculate markers and create plot with clustering_merged
  if("clustering_merged" %in% colnames(seurat@meta.data)) {
    
    Idents(seurat) <- "clustering_merged"
    
    #find markers per cluster
    markers_merged <- FindAllMarkers(seurat, logfc.threshold = 0.35, min.pct = 0.35, only.pos = T)
    
    top5_merged <- markers_merged %>%
      group_by(cluster) %>%
      top_n(n = 5, wt = avg_log2FC)
    top5_merged_unique <- unique(top5_merged$gene)
    
    results$markers_merged <- markers_merged
    
  }
  
  #proceed with clustering_fine calculations for all objects
  Idents(seurat) <- "clustering_fine"
  
  #find markers per cluster
  markers <- FindAllMarkers(seurat, logfc.threshold = 0.35, min.pct = 0.35, only.pos = T)
  
  top5 <- markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)
  top5_unique <- unique(top5$gene)
  
  results$markers <- markers

  #proceed with clustering_global calculations for all objects
  Idents(seurat) <- "clustering_global"
  
  #find markers per cluster
  markers <- FindAllMarkers(seurat, logfc.threshold = 0.35, min.pct = 0.35, only.pos = T)
  
  top5 <- markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)
  top5_unique <- unique(top5$gene)
  
  results$markers_global <- markers
  
  return(results)
  
}

#collect nr of genes from before and after filtering to monitor
gene_counts <- data.frame(Object = character(), Before_Filtering = integer(), After_Filtering = integer(), stringsAsFactors = FALSE)

for(i in files){
  obj_name <- strsplit(i, split = "_annotated_prefinal.RDS")[[1]]
  data <- readRDS(paste0(work_dir, '/', i))
  
  DefaultAssay(data) <- "RNA"
  
  #format metadata, assign age brackets depending on median
  data_formatted <- HAECA_metadata_formatting(data = data)
  
  print(table(data_formatted$Age, data_formatted$Age_brackets, useNA = 'always'))
  
  #collect for overview
  before_filtering <- nrow(data_formatted)
  
  #call tables after formatting
  print(paste0(obj_name, " after metadata formatting"))
  print(obj_name)
  print(colnames(data_formatted@meta.data))
  print(table(data_formatted$Age_numeric, useNA = 'always'))
  print(table(data_formatted$Age_brackets, useNA = 'always'))
  table(data_formatted$Decade, useNA = 'always')
  
#run filtering functions
data_filtered <- HAECA_gene_filtering(data = data_formatted, annotations = annotations, ambient_genes = ambient_genes)

#print results of filtered seurat object
print(obj_name)
print(colnames(data_filtered@meta.data))
print(data_filtered)
print(paste0("there are ", length(rownames(data_filtered)), " genes left in ", obj_name))

#collect for overview
after_filtering <- nrow(data_filtered)
gene_counts <- rbind(gene_counts, data.frame(Tissue = obj_name, before_filtering = before_filtering, after_filtering = after_filtering))

#save into output_dir
saveRDS(data_filtered, file = paste0(output_dir, "/", obj_name, "_filtered.RDS"))

#calculate top5 markers for global, merged and fine cluster
#call dotplot function
results <- HAECA_top5_markers(seurat = data)

#store outputs of function
#for merged
if("markers_merged" %in% names(results)) {
  markers_merged <- results$
  
  #save markers into output_dir
  write.csv(markers_merged, paste0(output_dir, "/", obj_name, "_clustering_merged_markers_filtered.csv"))
}

#for fine
markers <- results$markers

#save markers into output_dir
write.csv(markers, paste0(output_dir, "/", obj_name, "_clustering_fine_markers_filtered.csv"))

#for global
markers_global <- results$markers_global

#save markers into output_dir
write.csv(markers_global, paste0(output_dir, "/", obj_name, "_clustering_global_markers_filtered.csv"))

#prepare UMAPs of each filtered tissue object
#assign respective colors for clustering_fine
cols_fine <- tissue_color_mappings[[obj_name]]

#use clustering_merged as basis if present, otherwise proceed with clustering_fine
if("clustering_merged" %in% colnames(seurat@meta.data)) {
  seurat$clustering_fine <- seurat$clustering_merged
}

#Figure S1d - Export UMAPs
#Figure S4b - UMAP for liver with fine clusters used.
UMAP_clustering_fine <- plot(DimPlot(seurat, reduction = "umap", label = F, pt.size = 1, raster = F,
                                     cols = tissue_colors, group.by = "clustering_fine"))

UMAP_clustering_global <- plot(DimPlot(seurat, reduction = "umap", label = F, pt.size = 1, raster = F,
                                       cols = cols_global, group.by = "clustering_global"))

}

#save all before and after gene counts
write.csv(gene_counts, file = paste0(output_dir, "/gene_counts_filtering_summary.csv"), row.names = FALSE)S