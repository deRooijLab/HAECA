#Preprocessing of HAECA_44 - selection of ECs from single cell data
#source publication: Wu et al., 2024, 10.1038/s43587-024-00607-1
#downloaded from GEO under GSE255690

#open libraries
library(Seurat)
library(dplyr)
library(data.table)
library(tidyverse)
library(ggplot2)

#set working and output directory
work_dir <- setwd("~/HAECA/preprocessing/HAECA_44")
output_dir <- paste0(work_dir, "/output")
if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}

#open datafiles, add metadata, and generate Seurat object
data_dir <- paste0(work_dir, "/raw") #saved from GEO
files <- list.files(data_dir, full.names = TRUE, recursive = TRUE) 
files

#prepared metadata from GEO annotation and from supplementary tables from publication
metadata <- read.csv("HAECA_44_metadata.csv")

#make a vector with all sample names (= GEO accession numbers)
samples <- metadata$GEO_accession

#create subfolders with mtx, barcode & gene files for every sample
for(i in 1:length(samples)) {
  message(paste0('Processing ',samples[i]))
  patient <- list.files(data_dir, pattern = samples[i])
  patient_dir <- paste0(work_dir, "/aggregated/",samples[i])
  if(!dir.exists(patient_dir)){dir.create(patient_dir,recursive = T)}
  
  # custom function
  my_function <- function(x){
    file.rename( from = file.path(data_dir, x) ,
                 to = file.path(patient_dir, x) )
  }
  
  # apply the function to all files
  lapply(patient, my_function)
}

#process files, turn into individual Seurat objects, and add sample-level metadata
data_dir <- paste0(work_dir, "/aggregated/")
files <- list.files(data_dir, full.names = TRUE, recursive = FALSE)
files

#re-name the files in every subfolder to "matrix.mtx.gz", "barcodes.tsv.gz" and "features.tsv.gz"
for(i in 1:length(files)) {
  message(paste0('Processing ',files[i]))
  sample <- basename(files[i])
  file_dir <- file.path(data_dir,dir(data_dir)[grep(sample,dir(data_dir))])
  setwd(file_dir)
  matrix <- list.files(file_dir,pattern = '.mtx')
  new_name <- "matrix.mtx.gz"
  file.rename(from=matrix, to=new_name)
  genes <- list.files(file_dir,pattern = 'features')
  new_name <- "features.tsv.gz"
  file.rename(from=genes, to=new_name)
  barcodes <- list.files(file_dir,pattern = 'barcodes')
  new_name <- "barcodes.tsv.gz"
  file.rename(from=barcodes, to=new_name)
}

for(i in 1:length(files)) {
  message(paste0('Processing ',files[i]))
  sample <- basename(files[i])
  
  data <- Read10X(files[i])
  seurat <- CreateSeuratObject(data, project = sample, min.features = 100, min.cells = 3)
  
  sample_metadata <- metadata[which(metadata$GEO_accession == sample),]
  seurat$HAECA_ID <- sample_metadata$HAECA_ID
  seurat$Sample <- sample_metadata$Sample
  seurat$Donor <- sample_metadata$Donor
  seurat$Cell_source <- sample_metadata$Cell_source
  seurat$GEO_accession <- sample_metadata$GEO_accession
  seurat$Tissue <- sample_metadata$Tissue
  seurat$Age <- sample_metadata$Age
  seurat$Gender <- sample_metadata$Gender
  seurat$Ethnicity <- sample_metadata$Ethnicity
  seurat$BMI <- sample_metadata$BMI
  seurat$Tissue <- sample_metadata$Tissue
  seurat$Race <- sample_metadata$Race
  seurat$Smoking <- sample_metadata$Smoking
  seurat$Sample_info <- sample_metadata$Sample_info
  seurat$Method <- sample_metadata$Method
  seurat$Genome_build <- sample_metadata$Genome_build
  seurat$Sampling <- sample_metadata$Sampling
  
  save(seurat, file = paste0(output_dir, "/", sample,".RData"))
}

#combine all processed samples into one Seurat object, include joining layers
files <- list.files(output_dir, full.names = TRUE, recursive = TRUE) 
files 

for(i in 1:length(files)) {
  message(paste0('Processing',files[i]))
  
  # combine datasets
  if(i == 1)
  {
    get(load(files[i]))
    md <- seurat@meta.data
    sample_titles <- md$GEO_accession
    rownames(md) <- paste0(rownames(md), "_", sample_titles)
    seurat <- RenameCells(seurat, new.names = rownames(md))
    combined <- seurat
  }
  else 
  {
    get(load(files[i]))
    md <- seurat@meta.data
    sample_titles <- md$GEO_accession
    rownames(md) <- paste0(rownames(md), "_", sample_titles)
    seurat <- RenameCells(seurat, new.names = rownames(md))
    combined <- merge(combined, seurat)
    Layers(combined[["RNA"]])
    combined <- JoinLayers(combined)
  }}

combined #26174 features x 95032 samples
seurat <- combined

#quality control - check for number of genes per cell, and fraction of mito genes
Idents(seurat) <- "Donor"
seurat$percent.mt <- PercentageFeatureSet(object = seurat, pattern = "^MT-")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)

# Filter out low quality cells using selected thresholds, adapt for each experiment
#here percent.mt set to < 15 %
seurat <- subset(x = seurat, 
                 subset= 
                   (nFeature_RNA >= 200) & (nFeature_RNA <= 5000) &
                   (percent.mt < 15))

seurat #26174 features x 77474 samples

seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, selection.method = 'vst', nfeatures = 2000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, npcs = 50, verbose = FALSE)
ElbowPlot(seurat, ndims = 50)
DimHeatmap(seurat, dims = 1:12, cells = 500, balanced = TRUE)
DimHeatmap(seurat, dims = 13:24, cells = 500, balanced = TRUE)
DimHeatmap(seurat, dims = 25:36, cells = 500, balanced = TRUE)
#use 29 dims for clustering

seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:29, k.param = 25)
seurat <- FindClusters(seurat, resolution = 1.5)
seurat <- RunUMAP(seurat, dims = 1:29, min.dist = 0.3, n.neighbors = 35)
DimPlot(seurat, pt.size = 0.5, label = TRUE)
DimPlot(seurat, pt.size = 1, label = TRUE, group.by = "Age")
DimPlot(seurat, pt.size = 1, label = TRUE, group.by = "Donor")
DimPlot(seurat, pt.size = 1, label = TRUE, group.by = "Sample_info")

#prepare dotplot with celltype markers to pre-select ECs
DotPlot(seurat, features = c("PECAM1", "CDH5", "FLT1", "VWF", "KDR", "RGCC",
                             "CLDN5", "MMRN1", "PROX1", "PDPN", "LYVE1",
                             "PTPRC", "CD8A", "MRC1", "C1QA", 
                             "EPCAM", "DCN", "LUM", "RGS5",
                             "PDGFRB", "ACTA2", "TAGLN", "COL1A1",
                             "PROM1", "KIT", "CD34", "MSI2", "CD38"),
        cols = "RdBu") + RotatedAxis()

#EC markers
FeaturePlot(seurat, features = c("PECAM1","CDH5", "CLDN5", "KDR", "CD34", "PDGFRB"), cols = c('lightblue',"orange"))

#Lymphatic EC markers
FeaturePlot(seurat, features = c("PROX1", "PDPN", "LYVE1"), cols = c('lightblue',"orange"))

## set idents of EC clusters
EC <- subset(seurat, idents = c("11", "12", "13", "17", "22", "35")) 
#9898 ECs

#save subset
saveRDS(EC, file = paste0(output_dir, "/HAECA_44_selected_ECs.RDS"))
saveRDS(seurat, file = paste0(output_dir, "/HAECA_44_all_cells.RDS"))

sessionInfo()
