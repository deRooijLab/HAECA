**HAECA — Human Aging Endothelial Cell Atlas Analysis Pipelines**

This repository contains R scripts and workflows for processing and analyzing Human Aging Endothelial Cell Atlas (HAECA) single-cell RNA-seq data. The pipeline converts tissue-specific endothelial Seurat objects into pseudobulk datasets, and performs differential expression analysis. See bottom of this page for analyses on single-cell level (e.g., CDKN1A+ vs. CDKN1A- cells).

**Overview**

**Step 1: Create Pseudobulk Datasets**

Script: HAECA_pseudobulk_from_seurat

Input: Endothelial cell–filtered Seurat .RDS objects for each HAECA tissue

Data source: HAECA endothelial cell datasets (Zenodo: https://zenodo.org/records/16779452 - released upon publication)

What happens:

- Loading of all *_ECs_filtered.RDS files from work_dir.
- Aggregation of counts into pseudobulk datasets per Donor, Age_bracket, and HAECA_ID (categorized by single-cell and single-nucleus derived samples)
- Filtering out donor groups with fewer than 10 cells.
- Saving of pseudobulk .RDS files for downstream DESeq2 analysis.

Output bulk/bulk__age_brackets.rds


**Step 2: Wald Test (Differential Gene Expression Analysis)****

Script: HAECA_DGEA_aged_vs_young

Input: Pseudobulk .RDS files created in Step 1.

What happens:

- Loading of pseudobulk files for all tissues in work_dir.
- Creation of a DESeq2 dataset for each tissue.
- Wald tests (comparing aged vs. young), model: ~ Sampling_Study + Age_brackets

Output output/DGEA__wald.csv

**Installation**

Install dependencies:
install.packages(c("Seurat", "dplyr", "ggplot2", "stringr", "tidyr", "data.table", "tidyverse", "cowplot")) BiocManager::install(c("edgeR", "DESeq2", "EnhancedVolcano"))
