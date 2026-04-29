#HAECA - integration and annotation of ECs from Ovarium tissue
#Combine pre-selected ECs from HAECA_24, HAECA_44 and HAECA_45

#Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(data.table)
library(tidyverse)
library(viridis)
library(presto)
library(scDblFinder)

#####set directories#####
work_dir <- file.path("~/HAECA/tissue_analysis/Ovarium")
setwd(work_dir)
output_dir <- file.path("~/HAECA/global_analysis/tissues")
if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}

# Load data
sample_files <- list.files(work_dir, pattern = ".RDS", full.names = TRUE)
samples <- lapply(sample_files, readRDS)

#Check for specific patterns in gene names
gene_patterns <- c("ENS", "DARC")
pattern_check <- lapply(samples, function(sample) sapply(gene_patterns, function(pattern) any(str_detect(rownames(sample), pattern))))
print(pattern_check) # to verify presence
# ENS  DARC
# TRUE FALSE

# Merge Seurat objects iteratively
seurat <- merge(x = samples[[1]], y = samples[2:length(samples)])

#adjust metadata to uniform formatting
seurat$Tissue <- "Ovarium"
seurat$Cell_source <- "Ovarium"
seurat@meta.data <- seurat@meta.data %>%
  mutate(across(c(Method, Gender, Tissue, Sampling, Sample_info), str_to_title))

#combine layers
Layers(seurat[["RNA"]])
seurat <- JoinLayers(seurat)
seurat
#31708 features across 13593 samples within 1 assay

#repeat quality control
seurat$percent.mt <- PercentageFeatureSet(object = seurat, pattern = "^MT-")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
seurat <- subset(x = seurat, 
                 subset= 
                   (nFeature_RNA >= 200) & (nFeature_RNA <= 5000) &
                   (percent.mt < 15))

seurat #total 13562 ECs

#check number of ECs per donor, exclude donors < 100 ECs
metadata <- seurat@meta.data
remove <- metadata %>%
  dplyr::count(Donor) %>%
  filter(n < 100) %>%
  pull(Donor)
#Donor 3 has 81 ECs and will be removed 

Idents(seurat) <- "Donor"
seurat <- subset(x = seurat, idents = remove, invert = TRUE)
seurat@meta.data <- droplevels(seurat@meta.data)
seurat #total 13481 ECs

#### Normalize, scale, and run PCA ####
seurat <- seurat %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = 'vst', nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = TRUE)

ElbowPlot(seurat, ndims = 50)
DimHeatmap(seurat, dims = 1:12, cells = 500, balanced = TRUE)
DimHeatmap(seurat, dims = 13:24, cells = 500, balanced = TRUE)
DimHeatmap(seurat, dims = 25:36, cells = 500, balanced = TRUE)

# Find neighbors, clusters, and run UMAP
seurat <- seurat %>%
  FindNeighbors(dims = 1:24, k.param = 25) %>%
  FindClusters(resolution = 1,  cluster.name = "unintegrated_clusters") %>%
  RunUMAP(dims = 1:24, min.dist = 0.3, n.neighbors = 35, reduction.name = "umap.unintegrated")

# Generate DimPlots
dimplot_features <- c("unintegrated_clusters", "HAECA_ID", "Age", "Tissue", "Genome_build", "Sample_info", "Sampling")
DimPlot(seurat, reduction = "umap.unintegrated", group.by = dimplot_features, label = T)

# Generate DotPlot
marker_genes <- c("PECAM1", "CDH5", "FLT1", "VWF", "KDR", "RGCC", "CLDN5", "MMRN1", "PROX1", "PDPN", "PTPRC", "CD8A",
                  "MRC1", "C1QA", "CDH1", "EPCAM", "DCN", "LUM", "RGS5", "PDGFRB", "ACTA2", "TAGLN", "COL1A1",
                  "PROM1", "KIT", "CD34", "MSI2", "CD38", "TOP2A", "MKI67")
DotPlot(seurat, features = marker_genes, cols = "RdBu") + RotatedAxis()

#### Integration process ####
list <- SplitObject(seurat, split.by = "Donor")

# Normalize, find variable features, scale data, run PCA, and run scDblFinder for each dataset
list <- lapply(list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, verbose = TRUE)
  x <- RunPCA(x, verbose = TRUE)
  x <- as.SingleCellExperiment(x)
  x <- scDblFinder(x)
  x <- as.Seurat(x, assay = NULL)
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
                                  k.filter = 100, k.anchor = 5, reduction = "rpca")

seurat.integrated <- IntegrateData(anchorset = anchors, k.weight = 72, verbose = TRUE)
#k.weight adjusted

# Display the integrated object
seurat.integrated
#33708 features across 13481 samples within 2 assays 

#Note on scDoubletfinder: We recommend interpreting these scores in the context of cluster annotation. 
#All cells from a cluster with a large average doublet score should be considered suspect, and close neighbors of 
#problematic clusters should be treated with caution. A cluster containing only a small proportion of high-scoring 
#cells is safer, though this prognosis comes with the caveat that true doublets often lie immediately adjacent to their 
#source populations and end up being assigned to the same cluster. It is worth confirming that any interesting results of 
#downstream analyses are not being driven by those cells, e.g., by checking that DE in an interesting gene is not driven 
#solely by cells with high doublet scores. 

#run the dimensionality reduction
DefaultAssay(seurat.integrated) <- "integrated"

#standard workflow for visualization and clustering
seurat.integrated <- seurat.integrated %>%
  ScaleData(verbose = F) %>%
  RunPCA(npcs = 30, verbose = T)
ElbowPlot(seurat.integrated, ndims = 30)
DimHeatmap(seurat.integrated, dims = 13:24, cells = 500, balanced = TRUE)

# Find neighbors, clusters, and run UMAP
seurat.integrated <- seurat.integrated %>%
  RunUMAP(reduction = "pca", dims = 1:20,  reduction.name = "umap.integrated") %>%
  FindNeighbors(reduction = "pca", dims = 1:20, k.param = 25) %>%
  FindClusters(resolution = 1.5,  cluster.name = "integrated_clusters")

#check amount of predicted doublets
table(seurat.integrated$scDblFinder.class)
# doublet singlet 
# 617   12864

# Generate DimPlots
DimPlot(seurat.integrated, reduction = "umap", group.by = "integrated_clusters", label = T, pt.size = 1)
DimPlot(seurat.integrated, reduction = "umap", group.by = "scDblFinder.class",
        split.by = "scDblFinder.class")

#save QC plots
pdf(file = paste0(output_dir, "/UMAPs_post-integration.pdf"), width = 20, height = 12)
plot(DimPlot(seurat.integrated, group.by = c("integrated_clusters", "HAECA_ID", 
                                             "Age", "Donor",
                                             "Sampling", "scDblFinder.class"), pt.size = 1, label = T))
dev.off()

#save pre- and post-integration seurat objects
saveRDS(seurat, file = paste0(output_dir,"/ovarium_ECs_unintegrated.RDS"))
saveRDS(seurat.integrated, file = paste0(output_dir,"/ovarium_ECs_integrated.RDS"))

#### EC annotation ####
DefaultAssay(seurat.integrated) <- "RNA"

#Identify strong marker genes in each cluster
markers <- FindAllMarkers(seurat.integrated, logfc.threshold = 0.35, min.pct = 0.15,
                          only.pos = T, max.cells.per.ident = 500)
write.csv(markers, paste0(output_dir, "/clustering_round1.csv"))

#plot top-5 in dotplot
top5 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
top5_unique <- unique(top5$gene)

#plot top-5 in dotplot
pdf(file = paste0(output_dir, "/dotplot_top5_markers_integrated.pdf"), width = 20, height = 12)
DotPlot(seurat.integrated, features = top5_unique, cols = "RdBu", dot.scale = 5) + RotatedAxis() + 
  scale_colour_viridis(option = "D", discrete = FALSE) +
  theme(axis.text = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12))
dev.off()

#save QC plots
pdf(file = paste0(output_dir, "/UMAPs_EC_features_post-integration.pdf"), width = 20, height = 12)
plot(FeaturePlot(seurat.integrated, features = c("PECAM1", "CDH5", "ARL15", "HEY1","KDR", "RGCC", 
                                                 "ACKR1", "VCAM1","PROX1", "PDPN", "PLVAP", "ESM1"),
                 label = FALSE, pt.size = 1, cols = c('lightblue',"orange")))
dev.off()

#plot non-EC markers
pdf(file = paste0(output_dir, "/UMAPs_nonEC_features_post-integration.pdf"), width = 20, height = 12)
plot(FeaturePlot(seurat.integrated, features = c("PTPRC", "CD63", "DCN", "LUM", "MS4A1", "MS4A2", "EPCAM",
                                                 "CDH1", "PDGFRB", "ST18", "GFAP"),
                 label = FALSE, pt.size = 1, cols = c('lightblue',"orange")))
dev.off()

#plot QC scores
pdf(file = paste0(output_dir, "/UMAPs_QC_post-integration.pdf"), width = 20, height = 12)
plot(FeaturePlot(seurat.integrated, features = c("nFeature_RNA", "percent.mt","scDblFinder.score"),
                 label = FALSE, pt.size = 1, cols = c('lightblue',"orange")))
dev.off()

#plot vascular bed markers for annotation
DotPlot(seurat.integrated, features = c("GJA4", "GJA5", "FBLN2", "FBLN5", "HEY1", "DKK2", "CXCL12", "GLUL",
                                        "RGCC", "KDR", "FLT1", "CA4", "BTNL9",
                                        "VCAM1", "SELE", "SELP", "ACKR1", "NR2F2", "HDAC9",
                                        "ANGPT2", "ESM1", "PGF", "APLN", "LXN", "INSR", 
                                        "PROX1", "PDPN", "MMRN1", "CCL21", "TWIST2", "PTPRC", "PDGFRB")) + RotatedAxis() +
  scale_colour_viridis(option = "D", discrete = FALSE) +
  theme(axis.text = element_text(size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 11))

#check for common EC markers, features/cell and predicted doublets in each cluster
VlnPlot(seurat.integrated, features = "nFeature_RNA", pt.size = 0)
VlnPlot(seurat.integrated, features = "scDblFinder.score", pt.size = 0.5) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black")

#check how many doublets are classified per cluster
metadata <- seurat.integrated@meta.data
doublets <- metadata %>%
  filter(scDblFinder.class == "doublet") %>%
  dplyr::count(integrated_clusters, sort = TRUE)
#clusters 1, 0, 4 have most doublets

#find out frequency of doublets in each cluster
doublet_frequency <- metadata %>%
  group_by(integrated_clusters) %>%
  summarise(
    total = n(),
    doublet_count = sum(scDblFinder.class == "doublet")
  ) %>%
  mutate(doublet_frequency = doublet_count / total) %>%
  arrange(desc(doublet_frequency))
#clusters 22, 15, 13, 1 have > 7% doublets --> all will be removed.

#check markers for one specific cluster (if necessary)
markers11 <- FindMarkers(seurat.integrated, ident.1 = "11", only.pos = T)

#annotate the clusters to their predicted vascular bed subtype
Idents(seurat.integrated) <- "integrated_clusters"
levels(seurat.integrated)
old.ids <- c("0",
             "1", "2", "3", "4", "5",
             "6", "7", "8", "9", "10",
             "11", "12", "13", "14", "15",
             "16", "17", "18", "19", "20",
             "21", "22")

new.ids <- c("capillary_venous_1", 
             "stromal_remove",
             "vein_1",
             "low_quality_remove",
             "artery_1",
             "capillary_1", 
             "capillary_arterial",
             "artery_2",
             "vein_2",
             "lymphatic",
             "capillary_venous_2", 
             "low_quality_remove",
             "immune_cells_remove",
             "stromal_remove",
             "capillary_angiogenic",
             "stromal_remove",
             "stromal_remove",
             "vein_3",
             "low_quality_remove",
             "capillary_2", 
             "stromal_remove",
             "low_quality_remove",
             "doublets_remove"
)

names(new.ids) <- levels(seurat.integrated)
seurat.integrated <- RenameIdents(seurat.integrated, new.ids)
seurat.integrated <- AddMetaData(seurat.integrated, metadata = seurat.integrated@active.ident, col.name = "clustering_round1")
DimPlot(seurat.integrated, reduction = "umap", label = T, group.by = "clustering_round1")

#remove all non-EC clusters
Idents(seurat.integrated) <- "clustering_round1"
seurat.integrated_cleaned <- subset(seurat.integrated,
                                    idents = c("doublets_remove","immune_cells_remove",
                                               "stromal_remove","low_quality_remove"), invert = T)

seurat.integrated_cleaned #33708 features across 9089 samples within 2 assays

#order the clusters along the vascular bed, for nicer plotting
Idents(seurat.integrated_cleaned) <- "clustering_round1"
levels(seurat.integrated_cleaned)
DimPlot(seurat.integrated_cleaned, reduction = "umap", label = T, group.by = "clustering_round1")

#after removing clusters, re-run the scaling & dimensionality
DefaultAssay(seurat.integrated_cleaned) <- "integrated"

#standard workflow for visualization and clustering, and see if the integration worked well
seurat.integrated_cleaned <- seurat.integrated_cleaned %>%
  ScaleData(verbose = F) %>%
  RunPCA(npcs = 30, verbose = T)

ElbowPlot(seurat.integrated_cleaned, ndims = 30)
DimHeatmap(seurat.integrated_cleaned, dims = 1:12, cells = 500, balanced = TRUE)
DimHeatmap(seurat.integrated_cleaned, dims = 13:24, cells = 500, balanced = TRUE)

# Find neighbors, clusters, and run UMAP
seurat.integrated_cleaned <- seurat.integrated_cleaned %>%
  FindNeighbors(reduction = "pca", dims = 1:15, k.param = 25) %>%
  FindClusters(resolution = 1,  cluster.name = "integrated_clusters_cleaned") %>%
  RunUMAP(reduction = "pca", dims = 1:15, min.dist = 0.3, n.neighbors = 35)

# Generate DimPlots
DimPlot(seurat.integrated_cleaned, reduction = "umap", group.by = "integrated_clusters_cleaned", label = T, pt.size = 0.5)
DimPlot(seurat.integrated_cleaned, reduction = "umap", group.by = "clustering_round1", label = T, pt.size = 0.5)
DimPlot(seurat.integrated_cleaned, reduction = "umap", group.by = "scDblFinder.class")

#save QC plots
pdf(file = paste0(output_dir, "/UMAPs_integrated_2.pdf"), width = 20, height = 12)
plot(DimPlot(seurat.integrated_cleaned, group.by = c("integrated_clusters_cleaned",
                                                     "HAECA_ID", "Age", "Donor",
                                                     "Sampling", "scDblFinder.class"), pt.size = 1, label = T))
dev.off()

VlnPlot(seurat.integrated_cleaned, features = "scDblFinder.score", pt.size = 0.5)
VlnPlot(seurat.integrated_cleaned, features = "nFeature_RNA", pt.size = 0.5)

#repeat annotations after removing unwanted clusters
DefaultAssay(seurat.integrated_cleaned) <- "RNA"
markers <- FindAllMarkers(seurat.integrated_cleaned, logfc.threshold = 0.35, min.pct = 0.15,
                          only.pos = T, max.cells.per.ident = 500)
write.csv(markers, paste0(output_dir, "/clustering_round2.csv"))

top5 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
top5_unique <- unique(top5$gene)

DotPlot(seurat.integrated_cleaned, features = top5_unique, cols = "RdBu", dot.scale = 5) +
  RotatedAxis() + scale_colour_viridis(option = "D", discrete = FALSE) +
  theme(axis.text = element_text(size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 11))

DotPlot(seurat.integrated_cleaned, features = c("GJA4", "GJA5", "FBLN2", "FBLN5", "HEY1", "DKK2", "CXCL12", "GLUL",
                                                "RGCC", "KDR", "FLT1", "CA4", "BTNL9",
                                                "VCAM1", "SELE", "SELP", "ACKR1", "NR2F2", "HDAC9",
                                                "ANGPT2", "ESM1", "PGF", "APLN", "LXN", "INSR", 
                                                "PROX1", "PDPN", "MMRN1", "CCL21", "TWIST2", "PDGFRB")) + RotatedAxis() +
  scale_colour_viridis(option = "D", discrete = FALSE) +
  theme(axis.text = element_text(size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 11))

#save QC plots
pdf(file = paste0(output_dir, "/UMAPs_EC_features_integrated_2.pdf"), width = 20, height = 12)
plot(FeaturePlot(seurat.integrated_cleaned, features = c("PECAM1", "CDH5", "ARL15", "HEY1","KDR", "RGCC", 
                                                         "ACKR1", "VCAM1","PROX1", "PDPN", "PLVAP", "ESM1"),
                 label = FALSE, pt.size = 1, cols = c('lightblue',"orange")))
dev.off()

pdf(file = paste0(output_dir, "/UMAPs_nonEC_features_integrated_2.pdf"), width = 20, height = 12)
plot(FeaturePlot(seurat.integrated_cleaned, features = c("PTPRC", "CD63", "DCN", "LUM", "MS4A1", "MS4A2", "EPCAM",
                                                         "CDH1", "PDGFRB", "ST18", "GFAP"),
                 label = FALSE, pt.size = 1, cols = c('lightblue',"orange")))
dev.off()

old.ids <- c("0", "1", "2", "3", "4", "5",
             "6", "7", "8", "9", "10",
             "11")

new.ids <- c("capillary_venous_1", 
             "capillary_venous_2", 
             "artery_2",
             "vein_1",
             "capillary_arterial",
             "artery_1",
             "vein_2",
             "lymphatic",
             "capillary_1", 
             "capillary_angiogenic",
             "capillary_venous_3",
             "capillary_2")

names(new.ids) <- levels(seurat.integrated_cleaned)
seurat.integrated_cleaned <- RenameIdents(seurat.integrated_cleaned, new.ids)
seurat.integrated_cleaned <- AddMetaData(seurat.integrated_cleaned, metadata = seurat.integrated_cleaned@active.ident,
                                         col.name = "clustering_fine")

pdf(file = paste0(output_dir, "/UMAP_integrated_annotated_2.pdf"), width = 20, height = 12)
DimPlot(seurat.integrated_cleaned, reduction = "umap", group.by = "clustering_fine", label = T, raster = F)
dev.off()

#order the clusters along the vascular bed, going from artery - capillary - vein - lymphatic
Idents(seurat.integrated_cleaned) <- "clustering_fine"
levels(seurat.integrated_cleaned) <- c("artery_1",
                                       "artery_2",
                                       "capillary_arterial",
                                       "capillary_angiogenic",
                                       "capillary_1", 
                                       "capillary_2",
                                       "capillary_venous_1", 
                                       "capillary_venous_2", 
                                       "capillary_venous_3",
                                       "vein_1",
                                       "vein_2",
                                       "lymphatic")

#to change the visualization into the correct order, change clustering_fine to factor
seurat.integrated_cleaned$clustering_fine <- factor(seurat.integrated_cleaned$clustering_fine, levels = c("artery_1",
                                                                                                          "artery_2",
                                                                                                          "capillary_arterial",
                                                                                                          "capillary_angiogenic",
                                                                                                          "capillary_1", 
                                                                                                          "capillary_2",
                                                                                                          "capillary_venous_1", 
                                                                                                          "capillary_venous_2", 
                                                                                                          "capillary_venous_3",
                                                                                                          "vein_1",
                                                                                                          "vein_2",
                                                                                                          "lymphatic"))

cluster_order <- levels(seurat.integrated_cleaned)
pdf(file = paste0(output_dir, "/UMAP_post-integration_annotated.pdf"), width = 14, height = 12)
plot(DimPlot(seurat.integrated_cleaned, reduction = "umap", label = F, group.by = "clustering_fine", pt.size = 0.5,
             cols = c("#87092b", #artery_1
                      "#f70230", #artery_2
                      "#d18c9f", #capillary_arterial
                      "#D15C2C", #capillary_angiogenic
                      "#8cd1a9", #capillary_1
                      "#72a174", #capillary_2
                      "#16C4C1", #capillary_venous_1
                      "#6fc3ed", #capillary_venous_2
                      "#50c7c7",#capillary_venous_3
                      "#269feb", #vein_1
                      "#093c80",#vein_2
                      "#cfca42" #lymphatic
             )))
dev.off()

table(seurat.integrated_cleaned$clustering_fine)
#             artery_1             artery_2   capillary_arterial 
#                  846                  972                  888 
# capillary_angiogenic          capillary_1          capillary_2 
#                  372                  557                  155 
#   capillary_venous_1   capillary_venous_2   capillary_venous_3 
#                 1490                 1083                  259 
#               vein_1               vein_2            lymphatic 
#                  942                  825                  700

#export markers
markers <- FindAllMarkers(seurat.integrated_cleaned, logfc.threshold = 0.25, min.pct = 0.15, only.pos = T)
write.csv(markers, paste0(output_dir, "/clustering_final.csv"))

#calculate strong markers for dotplot
markers <- FindAllMarkers(seurat.integrated_cleaned, logfc.threshold = 0.35, min.pct = 0.35, only.pos = T)
top5 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

pdf(file = paste0(output_dir, "/Dotplot_post-integration_annotated.pdf"), width = 16, height = 5)
plot(DotPlot(seurat.integrated_cleaned, features = unique(top5$gene), cols = "RdBu") +
       RotatedAxis() + scale_colour_viridis(option = "D", discrete = FALSE))
dev.off()

#### merge clustering_fine to clustering_merged####
#Capillary_venous_2 + capillary_1
#vein1+2
levels(seurat.integrated_cleaned)
new.ids <- c("artery_1",
             "artery_2",
             "capillary_arterial",
             "capillary_angiogenic",
             "capillary_venous_2", 
             "capillary",
             "capillary_venous_1", 
             "capillary_venous_2", 
             "capillary_venous_3",
             "vein",
             "vein",
             "lymphatic")

names(new.ids) <- levels(seurat.integrated_cleaned)
seurat.integrated_cleaned <- RenameIdents(seurat.integrated_cleaned, new.ids)
seurat.integrated_cleaned <- AddMetaData(seurat.integrated_cleaned,
                                         metadata = seurat.integrated_cleaned@active.ident, col.name = "clustering_merged")

Idents(seurat.integrated_cleaned) <- "clustering_merged"
levels(seurat.integrated_cleaned) <- c("artery_1",
                                       "artery_2",
                                       "capillary_arterial",
                                       "capillary_angiogenic",
                                       "capillary",
                                       "capillary_venous_1", 
                                       "capillary_venous_2", 
                                       "capillary_venous_3",
                                       "vein",
                                       "lymphatic")

#to change the visualization into the correct order, change clustering_merged to factor
seurat.integrated_cleaned$clustering_merged <- factor(seurat.integrated_cleaned$clustering_merged,
                                                      levels = c("artery_1",
                                                                 "artery_2",
                                                                 "capillary_arterial",
                                                                 "capillary_angiogenic",
                                                                 "capillary",
                                                                 "capillary_venous_1", 
                                                                 "capillary_venous_2", 
                                                                 "capillary_venous_3",
                                                                 "vein",
                                                                 "lymphatic"))

pdf(file = paste0(output_dir, "/UMAP_post-integration_annotated_clustering_merged.pdf"), width = 14, height = 12)
plot(DimPlot(seurat.integrated_cleaned, reduction = "umap", label = F, group.by = "clustering_merged", pt.size = 0.5,
             cols = c("#87092b", #artery_1
                      "#f70230", #artery_2
                      "#d18c9f", #capillary_arterial
                      "#D15C2C", #capillary_angiogenic
                      "#8cd1a9", #capillary
                      "#16C4C1", #capillary_venous_1
                      "#6fc3ed", #capillary_venous_2
                      "#50c7c7",#capillary_venous_3
                      "#269feb", #vein
                      "#cfca42" #lymphatic
             )))
dev.off()

#calculate  markers from merged clustering
Idents(seurat.integrated_cleaned) <- "clustering_merged"
cluster_order <- levels(seurat.integrated_cleaned)

#set as factor to keep in row
seurat.integrated_cleaned$clustering_merged <- factor(seurat.integrated_cleaned$clustering_merged, levels = cluster_order)

markers <- FindAllMarkers(seurat.integrated_cleaned, logfc.threshold = 0.35, min.pct = 0.35, only.pos = T)
top5 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  as.data.frame()

pdf(file = paste0(output_dir, "/Dotplot_post-integration_annotated_merged.pdf"), width = 16, height = 5)
plot(DotPlot(seurat.integrated_cleaned, features = unique(top5$gene), cols = "RdBu") +
       RotatedAxis() + scale_colour_viridis(option = "D", discrete = FALSE))
dev.off()

#### global EC annotation ####
#this will allow a more harmonized analysis across tissues later on in the project
levels(seurat.integrated_cleaned)
new.ids <- c("artery",
             "artery",
             "capillary_arterial",
             "capillary_angiogenic",
             "capillary",
             "capillary_venous", 
             "capillary_venous", 
             "capillary_venous",
             "vein",
             "lymphatic")

names(new.ids) <- levels(seurat.integrated_cleaned)
seurat.integrated_cleaned <- RenameIdents(seurat.integrated_cleaned, new.ids)
seurat.integrated_cleaned <- AddMetaData(seurat.integrated_cleaned,
                                         metadata = seurat.integrated_cleaned@active.ident, col.name = "clustering_global")

pdf(file = paste0(output_dir, "/UMAP_post-integration_annotated_global.pdf"), width = 14, height = 12)
plot(DimPlot(seurat.integrated_cleaned, reduction = "umap", label = F, pt.size = 0.5, group.by = "clustering_global",
             cols = c("#87092b",
                      "#d18c9f",
                      "#D15C2C",
                      "#8cd1a9",
                      "#16C4C1",
                      "#269feb",
                      "#cfca42")))
dev.off()

#save seurat object, continue with visualization and analysis script
saveRDS(seurat.integrated_cleaned, paste0(output_dir, "/ovarium_ECs_annotated_prefinal.RDS"))

sessionInfo()