#global atlas integration by Sampling for Figures 1b, 1d and Figure S3a
#and by Tissue for Figure S1f (see below)

library(Seurat)
library(dplyr)
library(stringr)
library(tidyr)
library(data.table)
library(tidyverse)
library(stringr)
library(viridis)

#####set directories#####
work_dir <- file.path("~/HAECA/global_analysis/output")
setwd(work_dir)
output_dir <- file.path(paste0(work_dir, '/integrated'))
if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
setwd(output_dir)

#define colors for UMAPs
cols_global <- c("artery" = "#87092b",
                 "capillary_arterial" = "#d18c9f",
                 "capillary_angiogenic" = "#D15C2C",
                 "capillary" = "#8cd1a9",
                 "capillary_venous" = "#16C4C1",
                 "vein" = "#269feb",
                 "lymphatic" = "#cfca42")

tissue_col <- c(
  'Uterus' = "#d18c9f",
  'Brain' = "#4FAE62",
  'Lung' = "#F6C54D",
  'Breast' = "#E37D46",
  'Heart' = "#C02D45",
  'Muscle' = "#8ecae6",
  'Ovarium' = "#219ebc",
  'Liver' = "#023047",
  'Kidney' = "#dca100",
  'Skin' = "#482374FF",
  'Adipose' = "#fb8500",
  'Intestine' = "#24878EFF")

sampling_col <- c("#0D6153",
                  "#9F3864")

cols_brackets <- c(
  "young" = "#266c91",
  "aged" = "#f0cf16")

#load merged seurat object
seurat <- readRDS(paste0(work_dir, "/HAECA_global_pca.RDS"))

Idents(seurat) <- "Sampling"

# #### Integration process ####
list <- SplitObject(seurat, split.by = "Sampling")

#normalize, find variable features, scale data, run PCA
list <- lapply(list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, verbose = TRUE)
  x <- RunPCA(x, verbose = TRUE)
  x
})

# Select integration features and re-run scaling and PCA
features <- SelectIntegrationFeatures(object.list = list)

list <- lapply(list, function(x) {
  x <- ScaleData(x, features = features, verbose = TRUE)
  x <- RunPCA(x, features = features, verbose = TRUE)
  x
})

# Find anchors and integrate data
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features,
                                  k.filter = 200, k.anchor = 5, reduction = "rpca")

seurat.integrated <- IntegrateData(anchorset = anchors, k.weight = 200, verbose = TRUE)

# Display the integrated object
seurat.integrated

#run the dimensionality reduction after integration
DefaultAssay(seurat.integrated) <- "integrated"

#standard workflow for visualization and clustering
seurat.integrated <- seurat.integrated %>%
  ScaleData(verbose = F) %>%
  RunPCA(npcs = 30, verbose = T)

pdf(file = "integrated_elbowplot_heatmap.pdf", width = 14, height = 12)
ElbowPlot(seurat.integrated, ndims = 30)
DimHeatmap(seurat.integrated, dims = 1:12, cells = 500, balanced = TRUE)
DimHeatmap(seurat.integrated, dims = 13:24, cells = 500, balanced = TRUE)
dev.off()

# Find neighbors, clusters, and run UMAP
seurat.integrated <- seurat.integrated %>%
  RunUMAP(reduction = "pca", dims = 1:19,  reduction.name = "umap.integrated") %>%
  FindNeighbors(reduction = "pca", dims = 1:19, k.param = 25) %>%
  FindClusters(resolution = 1,  cluster.name = "integrated_clusters")

saveRDS(seurat.integrated, "global_integrated_sampling_pca.RDS")

#plot vascbed
Idents(seurat.integrated) <- "clustering_global"
levels(seurat.integrated)
seurat.integrated$clustering_global <- factor(seurat.integrated$clustering_global,levels = c("artery",
                                                                                             "capillary_arterial",
                                                                                             "capillary_angiogenic",
                                                                                             "capillary",
                                                                                             "capillary_venous",
                                                                                             "vein",
                                                                                             "lymphatic"))



pdf(file = paste0("Figure_1d_UMAP_vascbed_post-integration.pdf"), width = 14, height = 12)
plot(DimPlot(seurat.integrated, reduction = "umap.integrated", label = F, pt.size = 1, raster = F,
             group.by = "clustering_global", cols = cols_global))
dev.off()



pdf(file = "Figure_1b_UMAP_tissues_post-integration.pdf", width = 14, height = 12)
DimPlot(seurat.integrated, reduction = 'umap.integrated', pt.size = 1, label = F, raster = F, group.by = "Tissue",
        cols = tissue_col)

seurat.integrated$Age_brackets <- factor(seurat.integrated$Age_brackets,
                                         levels = c("young","aged"))

pdf(file = "Figure_S3a_UMAP_age-brackets_post-integration.pdf", width = 20, height = 10)
plot(DimPlot(seurat.integrated, reduction = 'umap.integrated', label = F, pt.size = 1, raster = F,
             cols = brackets_col, group.by = "Age_brackets", split.by = "Age_brackets"))
dev.off()

#global atlas integration by Tissue for Figure S1f
Idents(seurat) <- "Tissue"

# #### Integration process ####
list <- SplitObject(seurat, split.by = "Tissue")

# Normalize, find variable features, scale data, run PCA for each dataset
list <- lapply(list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, verbose = TRUE)
  x <- RunPCA(x, verbose = TRUE)
  x
})

# Select integration features and re-run scaling and PCA
features <- SelectIntegrationFeatures(object.list = list)

list <- lapply(list, function(x) {
  x <- ScaleData(x, features = features, verbose = TRUE)
  x <- RunPCA(x, features = features, verbose = TRUE)
  x
})

# Find anchors and integrate data
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features,
                                  k.filter = 200, k.anchor = 5, reduction = "rpca")

seurat.integrated <- IntegrateData(anchorset = anchors, k.weight = 200, verbose = TRUE)

#run the dimensionality reduction after integration
DefaultAssay(seurat.integrated) <- "integrated"

#standard workflow for visualization and clustering
seurat.integrated <- seurat.integrated %>%
  ScaleData(verbose = F) %>%
  RunPCA(npcs = 30, verbose = T)

saveRDS(seurat.integrated, "global_integrated_tissue.RDS")

pdf(file = "integrated_elbowplot_heatmap.pdf", width = 14, height = 12)
ElbowPlot(seurat.integrated, ndims = 30)
DimHeatmap(seurat.integrated, dims = 1:12, cells = 500, balanced = TRUE)
DimHeatmap(seurat.integrated, dims = 13:24, cells = 500, balanced = TRUE)
dev.off()

# Find neighbors, clusters, and run UMAP
seurat.integrated <- seurat.integrated %>%
  RunUMAP(reduction = "pca", dims = 1:20,  reduction.name = "umap.integrated") %>%
  FindNeighbors(reduction = "pca", dims = 1:20, k.param = 25) %>%
  FindClusters(resolution = 1,  cluster.name = "integrated_clusters")

saveRDS(seurat.integrated, "global_integrated_pca_tissue.RDS")

#plot vascbed
Idents(seurat.integrated) <- "clustering_global"
levels(seurat.integrated)
seurat.integrated$clustering_global <- factor(seurat.integrated$clustering_global,levels = c("artery",
                                                                        "capillary_arterial",
                                                                        "capillary_angiogenic",
                                                                        "capillary",
                                                                        "capillary_venous",
                                                                        "vein",
                                                                        "lymphatic"))

#Figure S1F
pdf(file = paste0("Figure_S1F_UMAP_vascbed_post-integration.pdf"), width = 14, height = 12)
plot(DimPlot(seurat.integrated, reduction = "umap.integrated", label = F, pt.size = 1, raster = F,
group.by = "clustering_global",
   cols = cols_global))
dev.off()

#FigureS1F-by-tissue
pdf(file = "Figure_S1F_UMAP_tissues_post-integration.pdf", width = 14, height = 12)
DimPlot(seurat.integrated, reduction = 'umap.integrated', pt.size = 1, label = F, raster = F, group.by = "Tissue",
        cols = tissue_col)
dev.off()