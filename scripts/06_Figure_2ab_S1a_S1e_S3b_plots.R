#Plot Figures 2a, 2b, S1a, S1e, S3b, data for supp. Table S3

library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(data.table)
library(tidyverse)
library(viridis)
library(presto)
library(enrichplot)
library(stringr)
library(tidyplots)
library(openxlsx)
library(ggExtra)
library(ggh4x)
library(ggpubr)

#####set directories#####
work_dir <- file.path("~/HAECA/global_analysis/output")
setwd(work_dir)
output_dir <- file.path(paste0(work_dir, '/plots'))
if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
setwd(output_dir)

# Define the colors for the bar graphs
cols <- c("#266c91","#f0cf16")

cols_brackets <- c(
  "young" = "#266c91",
  "aged" = "#f0cf16")

# Define the custom colors
cols_vascbed <- c(
  'artery' = "#87092b",
  'capillary_arterial' = "#d18c9f",
  'capillary_angiogenic' = "#D15C2C",
  'capillary' = "#8cd1a9",
  'capillary_venous' = "#16C4C1",
  'vein' = "#269feb",
  'lymphatic' = "#cfca42")

cols_tissue <- c(
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

#load merged seurat object
seurat <- readRDS("HAECA_global_pca.RDS")

#Bargraph Plots Figure S1a, Figure S3b
#prepare metadata
  metadata <- droplevels(seurat@meta.data)
  colnames(metadata)
  metadata$Age_brackets_donor <- paste0(metadata$Age_brackets, "_", metadata$Donor)
  metadata$Age_donor <- paste0(metadata$Age_numeric, "_", metadata$Donor)
  metadata$Age_brackets <- factor(metadata$Age_brackets, levels = c('young', 'aged'))
  
  #plot donors across tissues, colored by tissue, arranged by increasing nr
  FigS1a <- metadata %>%
    distinct(Donor, Tissue) %>%
    group_by(Tissue) %>%  
    summarize(Donor_Count = n(), .groups = "drop") %>%
    mutate(Tissue = fct_reorder(Tissue, Donor_Count, .desc = TRUE)) %>%
    tidyplot(x = Tissue, y = Donor_Count, color = Tissue) %>%
    add_barstack_absolute() %>%
    adjust_y_axis_title(title = "Number of Donors") %>%
    adjust_colors(new_colors = cols_tissue) %>%
    adjust_x_axis(rotate_labels = TRUE)

  #create bins for decades
  bins <- seq(10, 100, by = 10) - 1
  
  #plot histogram of number of donors across Age (numeric)
  FigS3b <- metadata %>%
    select("Age_donor", "Age_numeric", "Decade") %>%
    group_by(Age_donor) %>%
    distinct() %>%
    tidyplot(x = Age_numeric, color = Decade) %>%
    add_histogram(breaks = bins) %>%
    adjust_y_axis_title(title = "Number of Donors") %>%
    adjust_colors(new_colors = cols) %>%
    adjust_x_axis(breaks = bins)

#save plots to pdf
pdf(file = paste0(output_dir, "/Figures_S1a_S3b_bargraphs.pdf"), width = 10, height = 10)
plot(FigS1a)
plot(FigS3b)
dev.off()

#Table S3 - calculate age range and cutoff range for unique donors per tissue
  age_range_per_tissue <- metadata %>%
    group_by(Tissue, Age_brackets)%>%
    summarise(
      min_Age = min(Age_numeric, na.rm = TRUE),
      max_Age = max(Age_numeric, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(id = paste0(Tissue, "_", Age_brackets)) %>%
    select(-Tissue, -Age_brackets)
  
  metadata_unique <- metadata %>%
    select(Tissue, Age_brackets, Age_numeric, median_age, cutoff_age, Donor) %>%
    group_by(Tissue, Age_brackets) %>%
    distinct(Donor, .keep_all = T) %>%
    mutate(count = n()) %>%
    mutate(id = paste0(Tissue, "_", Age_brackets)) %>%
    left_join(y = age_range_per_tissue, by = 'id') %>%
    select(-id)
  
  #Supplementary Table S3
  write.csv(metadata_unique, paste0(output_dir, '/', 'Table_S3_Age_range_per_tissue.csv'))
  
#Fig2a - individual tissue color added to header strip
  Fig2a <- ggplot(metadata_unique,
                aes(x = factor(Age_brackets, levels = c("young", "aged")), 
                    y = Age_numeric, fill = Age_brackets)) +
    geom_boxplot() +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    geom_hline(data = metadata_unique, aes(yintercept = median_age), linetype = "dashed", color = "#4b4e52") +
    geom_hline(data = metadata_unique, aes(yintercept = cutoff_age), linetype = "dashed", color = "darkred") +
    geom_text(data = metadata_unique, aes(x = Age_brackets, y = 5, 
                                          label = count), color = "black", size = 4) +
    theme_bw() +
    labs(y = "Age (Years)", fill = "Age Bracket") +
    scale_fill_manual(values = cols_bracket) +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      legend.position = "none",
      axis.title.x = element_blank()
    ) +
    facet_grid(~ Tissue, scales = "free") +
    theme(
      strip.text = element_text(size = 12, face = "bold", color = "white"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      legend.position = "none",
      axis.title.x = element_blank()) +
    facet_grid2(~ Tissue, strip = strip_themed(
      background_x = list(element_rect(fill = "#fb8500"),
                          element_rect(fill = "#4FAE62"),
                          element_rect(fill = "#E37D46"),
                          element_rect(fill = "#C02D45"),
                          element_rect(fill = "#24878EFF"),
                          element_rect(fill = "#dca100"),
                          element_rect(fill = "#023047"),
                          element_rect(fill = "#F6C54D"),
                          element_rect(fill = "#8ecae6"),
                          element_rect(fill = "#219ebc"),
                          element_rect(fill = "#482374FF"),
                          element_rect(fill = "#d18c9f")
      )))

pdf(file = paste0(output_dir, "/Figure_2a_bargraph_age_range.pdf"), width = 20, height = 5)
plot(Fig2a)
dev.off()

#Figure S1e - Create contingency table for Age_brackets, Donor, and clustering_global
df_brackets <- table(seurat$Tissue, seurat$Age_brackets, seurat$Donor, seurat$clustering_global) %>%
as.data.frame() %>%
filter(Freq > 0)
  
  # Rename columns
  colnames(df_brackets) <- c("Tissue", "Age", "Donor", "clustering_global", "count")
  
  # Ensure clustering_global follows the correct order
  df_brackets$clustering_global <- factor(df_brackets$clustering_global, levels = names(cols_vascbed))
  
  # Calculate percentage per donor
  df_brackets_percentage <- df_brackets %>%
    group_by(Tissue, Donor) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup() %>%
    group_by(Tissue, clustering_global, Age) %>%
    mutate(
      median_percentage = median(percentage),
      mean_percentage = mean(percentage),
      sd_percentage = sd(percentage),
      sem_percentage = sd(percentage) / sqrt(n()))  %>%
    group_by(clustering_global, Age) %>%
    mutate(
      median_percentage_global = median(percentage),
      mean_percentage_global = mean(percentage),
      sd_percentage_global = sd(percentage),
      sem_percentage_global = sd(percentage) / sqrt(n())) %>%
    mutate(Age = factor(Age, levels = c("young", "aged")))
  
  #prepare workbook with all necessary dataframes
  data <- list(
    df_brackets = df_brackets, df_brackets_percentage = df_brackets_percentage,
    df_fold_change_tissue = df_fold_change_tissue
  )
  
  wb <- createWorkbook()
  
  for(i in names(data)) {
    addWorksheet(wb, i)
    writeData(wb, i, data[[i]])
  }

#Supplementary Tables SNN and SNN
saveWorkbook(wb, file = paste0(output_dir, "/HAECA_folchange_data.xlsx"), overwrite = TRUE)

#Figure S1e - plot bargraphs by tissue
#Create a contingency table of Tissue, HAECA_ID, and clustering_global
df_bargraphs <- table(seurat$Tissue, seurat$HAECA_ID, seurat$Donor, seurat$clustering_global) %>%
    as.data.frame()
  
  # Rename the columns
  colnames(df_bargraphs) <- c("Tissue", "HAECA_ID", "Donor", "clustering_global", "count")
  
  # Convert clustering_global to a factor with the specified levels
  df_bargraphs$clustering_global <- factor(df_bargraphs$clustering_global, levels = names(cols_vascbed))
  
  # Calculate percentages
  df_bargraphs_percentage <- df_bargraphs %>%
    filter(count > 0) %>%
    group_by(Tissue, HAECA_ID) %>%
    mutate(percentage_study = count / sum(count) * 100) %>%
    group_by(Tissue, Donor) %>%
    mutate(percentage_donor = count / sum(count) * 100)
  
  #Figure S1e - plot Percentages with Grouping by Tissue and Donor
  FigS1e <- ggplot(df_bargraphs_percentage, aes(x = Donor, y = percentage_donor, fill = clustering_global)) +
    geom_bar(stat = "identity") +
    facet_grid(Tissue ~ ., scales = "free_y", space = "free_y") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Percentage per Donor, grouped by Tissue",
         x = "Donor",
         y = "Percentage") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = cols_vascbed) +
    theme_classic() +
    theme(
      strip.text.y = element_text(angle = 0),
      panel.spacing.y = unit(0, "lines"),
      axis.text.x = element_text(angle = 90, hjust = 1))

pdf(file = paste0(output_dir, "/Figure_S1e_bargraph_by_donor.pdf"), width = 10, height = 30)
plot(FigS1e)
dev.off()



#Figure 2b - violin foldchange by cluster
# Create contingency table for Age_brackets, Donor, and clustering_global
df_brackets <- table(seurat$Tissue, seurat$Age_brackets, seurat$Donor, seurat$clustering_global) %>%
  as.data.frame() %>%
  filter(Freq > 0)

# Rename columns
colnames(df_brackets) <- c("Tissue", "Age", "Donor", "clustering_global", "count")

# Ensure clustering_global follows the correct order
df_brackets$clustering_global <- factor(df_brackets$clustering_global, levels = names(cols_vascbed))

# Calculate percentage per donor
df_brackets_percentage <- df_brackets %>%
  group_by(Tissue, Donor) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup() %>%
  group_by(Tissue, clustering_global, Age) %>%
  mutate(
    median_percentage = median(percentage),
    mean_percentage = mean(percentage),
    sd_percentage = sd(percentage),
    sem_percentage = sd(percentage) / sqrt(n()))  %>%
  group_by(clustering_global, Age) %>%
  mutate(
    median_percentage_global = median(percentage),
    mean_percentage_global = mean(percentage),
    sd_percentage_global = sd(percentage),
    sem_percentage_global = sd(percentage) / sqrt(n())) %>%
  mutate(Age = factor(Age, levels = c("young", "aged")))

# Calculate fold_change of median percentage across clusters and within tissue
df_fold_change_tissue <- df_brackets_percentage %>%
  select(Tissue, Age, clustering_global, median_percentage) %>%
  distinct() %>%
  pivot_wider(names_from = Age, values_from = median_percentage) %>%
  mutate(
    fold_change = aged / young,
    log2_fold_change = log2(fold_change)
  ) %>%
  group_by(clustering_global) %>%
  arrange(log2_fold_change) %>%
  mutate(rank_cluster = row_number()) %>%
  ungroup() %>%
  group_by(Tissue) %>%
  arrange(log2_fold_change) %>%
  mutate(rank_tissue = row_number()) %>%
  mutate(clustering_global = factor(clustering_global, levels = names(cols_vascbed)))

df_fold_change_tissue <- df_fold_change_tissue %>%
  mutate(clustering_global = factor(clustering_global, levels = names(cols_vascbed)))

p5 <- df_fold_change_tissue %>%
  tidyplot(x = log2_fold_change, y = clustering_global, color = clustering_global) %>% 
  add_data_points_jitter() %>%
  add_reference_lines(x = 0) %>% 
  add_violin() %>%
  adjust_x_axis_title(title = 'median_log2_fold_change') %>%
  adjust_colors(new_colors = cols_vascbed) %>%
  reorder_y_axis_labels(levels(df_fold_change_tissue$clustering_global)) %>%
  remove_legend() %>%
  remove_y_axis_title()

#Supplementary Tables SNN and SNN
#prepare workbook with all necessary dataframes
data <- list(
  df = df, df_percentage = df_percentage, df_donor = df_donor, 
  df_bargraphs = df_bargraphs, df_bargraphs_percentage = df_bargraphs_percentage, 
  df_clusters = df_clusters, df_clusters_percentage = df_clusters_percentage
)

wb <- createWorkbook()

for(i in names(data)) {
  addWorksheet(wb, i)
  writeData(wb, i, data[[i]])
}

saveWorkbook(wb, file = paste0(output_dir, "/HAECA_visualization_data.xlsx"), overwrite = TRUE)