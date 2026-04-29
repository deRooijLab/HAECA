#Merge and cluster global atlas from filtered tissues

library(Seurat)
library(dplyr)
library(stringr)
library(tidyr)
library(data.table)
library(tidyverse)
library(presto)
library(pheatmap)
library(viridis)

#####set directories#####
tissue_dir <- file.path("~/HAECA/global_analysis/tissues/filtered")
work_dir <- file.path("~/HAECA/global_analysis/output")
if(!dir.exists(work_dir)){dir.create(work_dir, recursive = T)}
setwd(work_dir)

#define color scheme for tissues
tissue_col <- c(
  'Adipose' = "#fb8500",
  'Brain' = "#4FAE62",
  'Breast' = "#E37D46",
  'Heart' = "#C02D45",
  'Intestine' = "#24878EFF",
  'Kidney' = "#dca100",
  'Liver' = "#023047",
  'Lung' = "#F6C54D",
  'Muscle' = "#8ecae6",
  'Ovarium' = "#219ebc",
  'Skin' = "#482374FF",
  'Uterus' = "#d18c9f")

#load individual filtered tissues
sample_files <- list.files(tissue_dir, pattern = ".RDS", full.names = TRUE)
samples <- lapply(sample_files, readRDS)
sample_names <- basename(sample_files)

#merge Seurat objects
seurat <- merge(x = samples[[1]], y = samples[2:length(samples)], project = "HAECA_global")

#combine layers
Layers(seurat[["RNA"]])
seurat <- JoinLayers(seurat)
seurat@meta.data <- droplevels(seurat@meta.data)

#report metadata contents
colnames(seurat@meta.data)

#export metadata for table S1
metadata <- seurat@meta.data
colnames(metadata)

meta <- metadata[,-c(1:3, 16:20, 22, 27)] %>%
  group_by(Tissue, HAECA_ID, Donor) %>%
  distinct()

write.csv(meta, "global_metadata_table_S1.csv", row.names = FALSE, quote = FALSE)

#normalize, scale, and run PCA
seurat <- seurat %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = 'vst', nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = TRUE)

pdf(file = "unintegrated_elbowplot_heatmaps.pdf", width = 14, height = 12)
plot(ElbowPlot(seurat, ndims = 50))
DimHeatmap(seurat, dims = 1:12, cells = 500, balanced = TRUE)
DimHeatmap(seurat, dims = 13:24, cells = 500, balanced = TRUE)
DimHeatmap(seurat, dims = 25:38, cells = 500, balanced = TRUE)
dev.off()

#find neighbors, clusters, and run UMAP
DefaultAssay <- "RNA"

seurat <- seurat %>%
  FindNeighbors(dims = 1:20, k.param = 25) %>%
  FindClusters(resolution = 1,  cluster.name = "unintegrated_clusters") %>%
  RunUMAP(dims = 1:20, min.dist = 0.3, n.neighbors = 35, reduction.name = "umap.unintegrated")

#save seurat object
saveRDS(seurat, paste0("HAECA_global_pca.RDS"))

#prepare UMAPs
pdf(file = "UMAPs_pre-integration.pdf", width = 26, height = 12)
plot(DimPlot(seurat, reduction = "umap.unintegrated", group.by = c("unintegrated_clusters","Tissue",
                                                                   "Sampling"), label = T, raster = T))
dev.off()

seurat$Tissue <- as.character(seurat$Tissue)
Idents(seurat) <- "Tissue"

#Figure S1C - export UMAP with Sampling
pdf(file = "UMAP_Sampling.pdf", width = 14, height = 12)
DimPlot(seurat, pt.size = 1, label = F, raster = F, group.by = "Sampling", 
        cols = c("#0D6153", "#9F3864"))
dev.off()

### marker calculation ####
DefaultAssay(seurat) <- "RNA"

#calculate  Tissue-specific markers
Idents(seurat) <- "Tissue"
markers <- FindAllMarkers(seurat, logfc.threshold = 0.35, min.pct = 0.35, only.pos = T)
top5 <- markers %>%
  group_by(cluster) %>%
  filter(!duplicated(gene)) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  as.data.frame()

write.csv(markers, "global_markers.csv")

#calculate  Cluster-specific markers
seurat$Tissue_cluster <- paste0(seurat$Tissue, "_", seurat$clustering_global)

Idents(seurat) <- "Tissue_cluster"
markers <- FindAllMarkers(seurat, logfc.threshold = 0.35, min.pct = 0.35, only.pos = T)

write.csv(markers, "tissue_cluster_markers.csv")

top3 <- markers %>%
  group_by(cluster) %>%
  filter(!duplicated(gene)) %>%
  top_n(n = 3, wt = avg_log2FC) %>%
  as.data.frame()

#change formatting
seurat$Tissue_cluster <- gsub("_", "-", seurat$Tissue_cluster)
expr_matrix <- AverageExpression(seurat, features = unique(top3$gene), group.by = "Tissue_cluster")$RNA

#scale the expression values
expr_matrix_scaled <- t(scale(t(expr_matrix)))

#add annotation
  dplyr::select(Tissue_cluster, Tissue, clustering_global) %>%
  distinct()

tissues_vector <- tissue_lookup$Tissue[match(colnames(expr_matrix_scaled), tissue_lookup$Tissue_cluster)]
clusters_vector <- tissue_lookup$clustering_global[match(colnames(expr_matrix_scaled), tissue_lookup$Tissue_cluster)]

annotation_col <- data.frame(
  Tissue = tissues_vector,
  Cluster = clusters_vector
)

rownames(annotation_col) <- colnames(expr_matrix_scaled)

annotation_colors <- list(
  Tissue = tissue_col,
  Cluster = cols_global
)

#Figure 1C - heatmap
pdf(file = paste0("Figure_1c_Heatmap_Tissue_cluster_markers_top3.pdf"), width = 10, height = 12)
pheatmap(expr_matrix_scaled,
         color = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         fontsize_row = 6,
         fontsize_col = 8,
         border_color = NA)
dev.off()