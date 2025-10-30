#script create bulk from each tissue
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(data.table)
library(tidyverse)

work_dir <- file.path("~/HAECA/global_analysis/tissues/filtered")
setwd(work_dir)
output_dir <- paste0("~/HAECA/global_analysis/bulk")
if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}

# Load data
files <- list.files(work_dir, pattern = "_ECs_filtered.RDS", full.names = FALSE)
tissues <- str_remove(files, "_ECs_filtered.RDS")

donor_counts <- data.frame() 
donor_cell_counts <- data.frame()

#loop through each tissue seurat object:
for(tissue in tissues) {
  
  file_name <- paste0(tissue, "_ECs_filtered.RDS")
  seurat <- readRDS(file_name)
  seurat$Sampling_Study <- paste0(seurat$Sampling, "_", seurat$HAECA_ID)
  
  #create bulk object split by Age brackets
  if (tissue == "kidney") {
    grouping_vars <- c("Donor", "Age_brackets") #only one study present in Kidney tissue
  } else {
    grouping_vars <- c("Donor", "Age_brackets", "Sampling_Study")}

  #prepare table for number of donors per age bracket and sampling_study
  donor_count <- seurat@meta.data %>%
    group_by(Age_brackets, Sampling_Study) %>%
    summarise(Donor_Count_per_bracket_study = n_distinct(Donor), .groups = "drop") %>%
    group_by(Age_brackets) %>%
    mutate(Donor_Count_per_bracket = sum(Donor_Count_per_bracket_study)) %>%
    ungroup() %>%
    mutate(Tissue = str_to_title(tissue))
  
  donor_counts <- bind_rows(donor_counts, donor_count)
  
  donor_cell_count <- seurat@meta.data %>%
    select(all_of(grouping_vars)) %>%
    mutate(Group = do.call(paste, c(., sep = "_"))) %>%
    count(Donor, Group, name = "CellCount") %>%
    mutate(Tissue = str_to_title(tissue))
  
  donor_cell_counts <- bind_rows(donor_cell_counts, donor_cell_count)
  
  #remove donor group with < 10 cells before bulk object
  filtered_groups <- donor_cell_count %>%
    filter(CellCount >= 10) %>%
    pull(Group)
  
  #change to grouping vars
  seurat$Group <- seurat@meta.data %>%
    dplyr::select(all_of(grouping_vars)) %>%
    mutate(Group = do.call(paste, c(., sep = "_"))) %>%
    pull(Group)
  
  #subset based on filtered donors containing >10 cells each
  Idents(seurat) <- "Group"
  seurat <- subset(seurat, idents = filtered_groups)

  bulk <- AggregateExpression(seurat, assays = "RNA", return.seurat = TRUE,
                      group.by = grouping_vars)
  
  #save bulk object
  output_file <- file.path(output_dir, paste0("bulk_", tissue, ".rds"))
  saveRDS(bulk, file = output_file)
  
}

write.csv(donor_counts, file = paste0(output_dir, "/donor_counts.csv"), row.names = F)
write.csv(donor_cell_counts, file = paste0(output_dir, "/donor_cell_counts.csv"), row.names = F)
