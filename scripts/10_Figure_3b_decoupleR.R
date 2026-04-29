#TF activity inference HAECA for Figure 3b

#### open libraries ####
library(Seurat)
library(decoupleR)
library(OmnipathR)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(qvalue)

#source directory of wald test results
source_dir <- paste0("~/HAECA/global_analysis/deseq2")
setwd(source_dir)
output_dir <- paste0("~/HAECA/global_analysis/decoupler")
if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}

#get function
#decoupler run
run_decoupler <- function(data_path, tissue, output_dir, net, n_tfs = 25) {
  
  # Load DGEA test results
  if (!file.exists(data_path)) {
    message("File not found: ", data_path)
    return(NULL)
  }
  
  DGEA <- read.csv(data_path, row.names = 1) %>%
    remove_rownames() %>%
    column_to_rownames(var = "gene")
  
  mat <- DGEA[, 'stat', drop = FALSE] %>%
    na.omit()
  
  # Run ULM
  contrast_acts <- run_ulm(mat = mat, net = net, .source = 'source', .target = 'target',
                           .mor = 'mor', minsize = 5)
  
  # Save contrast_acts results
  tfact_file <- file.path(output_dir, paste0("TFact_", tissue, ".csv"))
  write.csv(contrast_acts, file = tfact_file)

  return(contrast_acts)
  
}

# Load data
files <- list.files(pattern = "DGEA_", full.names = FALSE)
tissues <- unique(str_extract(files, "(?<=DGEA_)[a-zA-Z]+"))

#### CollecTRI network ####
omnipath_cache_wipe()
net <- decoupleR::get_collectri(organism = 'human', split_complexes = FALSE)

results_list <- list()

#loop through each tissue and condition, perform analysis on each combination and save pdf with plots
for (tissue in tissues) {
  
  tissue_dir <- file.path(output_dir, tissue)
  if (!dir.exists(tissue_dir)) dir.create(tissue_dir, recursive = TRUE)
  
  data_path <- file.path(source_dir, paste0("DGEA_", tissue, "_wald.csv"))
    
  #run decoupler analysis
  results <- run_decoupler(data_path, tissue, output_dir = tissue_dir, net, n_tfs = 30)
  
  #collect plots and results
  results_list[[tissue]] <- results$contrast_acts
   
}

####combined analysis across tissues #####


data <- bind_rows(results_list, .id = "Tissue")

#add adjusted pvalue to reults
df$p_adj <- p.adjust(df$p_value, method = "BH")

#save for Supplementary Table #9
tfact_file <- file.path(output_dir, paste0("/Table_S9_HAECA.csv"))
write.csv(df, file = tfact_file, row.names = F, quote = F)

#For Figure 3b - define TFs of interest
selected_tfs <- c("YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4")

df_filtered <- df %>%
  filter(source %in% selected_tfs) %>%
  mutate(
    Significant = p_value < 0.05,
    log_p = -log10(p_value)
  )

df_filtered <- df_filtered %>%
  mutate(score_capped = pmin(pmax(score, -3), 3))

df_filtered$log_p_scaled <- scales::rescale(df_filtered$log_p, to = c(0, 5))

# Final plot
pdf(file = "Figure_3b_dotplot_mechanoTFs.PDF", width = 4, height = 2)
ggplot(df_filtered, aes(x = Tissue, y = source)) +
  geom_point(
    aes(
      size = log_p,
      fill = score_capped,
      color = Significant,
      stroke = ifelse(Significant, 1, 0.2)
    ),
    shape = 21,
    alpha = 0.8
  ) +
  scale_size_continuous(range = c(0, 10), limits = c(0, 5)) +
  scale_fill_viridis(option = "viridis", name = "Score") +  # yellow high, blue low
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  theme_classic() +
  labs(
    title = "Mechanosensitive TFs",
    x = "Tissue",
    y = "Transcription Factor",
    size = "-log10(p)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )
dev.off()