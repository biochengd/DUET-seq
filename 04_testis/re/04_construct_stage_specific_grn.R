# ==============================================================================
# DUET-SEQ ANALYSIS PIPELINE
# Script: 04_construct_stage_specific_grn.R
# Description: Construction of stage-specific Gene Regulatory Networks (GRNs)
#              using FigR predictions and AUCell activity scores.
#              Compares Round Spermatids vs Elongating Spermatids.
# ==============================================================================

# 1. Setup and Configuration ---------------------------------------------------

# Load necessary libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(FigR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(igraph)
  library(ggraph)
  library(qs)
  library(glue)
  library(tibble)
})

# Define Paths (Anonymized)
DATA_DIR   <- "./data"       # Directory containing Seurat object and FigR results
OUTPUT_DIR <- "./output/GRN" # Directory for results
dir.create(file.path(OUTPUT_DIR, "plots"), showWarnings = FALSE, recursive = TRUE)

# Define Analysis Parameters
FIGR_SCORE_THRESHOLD <- 1.0   # Threshold for FigR regulation score
TF_ACTIVITY_THRESHOLD <- 0.05 # Minimum average TF activity to consider active
NETWORK_QUANTILE <- 0.5       # Keep top 50% links based on combined activity
MIN_TARGETS_FOR_PLOT <- 3     # Minimum targets a TF must have to appear in plots
TOP_N_EDGES_PLOT <- 150       # Number of edges to visualize

# 2. Data Loading and Preprocessing --------------------------------------------

message("Loading data...")

# Load Seurat object (RNA + AUCell)
# Expects an object with 'AUCell' assay
seu <- qread(file = glue("{DATA_DIR}/seurat_aucell_processed.qs"))

# Load FigR results
figR.d <- readRDS(file = glue("{DATA_DIR}/figR_results.rds"))

# Refine Cell Types (Merge subtypes for robust GRN construction)
# Logic: Meig1/Crem -> Early Elongating; Akap4/Odf1 -> Late Elongating
seu$celltype_refined <- as.character(seu$celltype_re_lvl2)

seu$celltype_refined[seu$celltype_refined %in% c("Elongating_Spermatids_Meig1", "Elongating_Spermatids_Crem")] <- "Elongating_Spermatids_Early"
seu$celltype_refined[seu$celltype_refined %in% c("Elongating_Spermatids_Akap4", "Elongating_Spermatids_Odf1")] <- "Elongating_Spermatids_Late"

# Set factor levels
refined_levels <- c("Round_Spermatids_Early", "Round_Spermatids_Late", 
                    "Elongating_Spermatids_Early", "Elongating_Spermatids_Late")
seu$celltype_refined <- factor(seu$celltype_refined, levels = refined_levels)

# 3. Core TF Filtering ---------------------------------------------------------

message("Filtering high-quality regulatory connections...")

# Get TFs available in AUCell matrix
rasMat <- GetAssayData(seu, assay = "AUCell", layer = "data")
tfs_in_aucell <- rownames(rasMat)

# Filter FigR data based on score and intersection with AUCell TFs
figR.d_filtered <- figR.d %>%
  dplyr::filter(Score > FIGR_SCORE_THRESHOLD) %>%
  dplyr::filter(Motif %in% tfs_in_aucell)

# Identify Core TFs
core_tfs <- unique(figR.d_filtered$Motif)
message(glue("Identified {length(core_tfs)} core TFs for network construction."))

# Save Core TF list
write.csv(core_tfs, file = glue("{OUTPUT_DIR}/core_tfs_list.csv"), row.names = FALSE)

# 4. TF Activity Calculation ---------------------------------------------------

message("Calculating stage-specific TF activity...")

# Extract metadata
cell_groups <- seu@meta.data %>%
  dplyr::select(celltype_refined) %>%
  tibble::rownames_to_column("cell")

# Calculate average TF activity per cell type
tf_activity_stats <- as.data.frame(t(rasMat)) %>%
  tibble::rownames_to_column("cell") %>%
  tidyr::pivot_longer(-cell, names_to = "TF", values_to = "activity") %>%
  dplyr::left_join(cell_groups, by = "cell") %>%
  dplyr::group_by(celltype_refined, TF) %>%
  dplyr::summarise(
    avg_activity = mean(activity),
    median_activity = median(activity),
    .groups = 'drop'
  )

# 5. Network Construction Function ---------------------------------------------

create_condition_network <- function(condition_name, 
                                     network_data, 
                                     activity_stats, 
                                     act_thresh = 0.05, 
                                     quantile_thresh = 0.5) {
  
  # Filter for TFs active in this specific condition
  condition_tf_activity <- activity_stats %>%
    filter(celltype_refined == condition_name) %>%
    rename(tf_avg_activity = avg_activity)
  
  active_tfs <- condition_tf_activity %>%
    filter(tf_avg_activity > act_thresh) %>%
    pull(TF)
  
  # Construct network
  # Link Activity = TF Activity * log2(FigR Score + 1)
  condition_network <- network_data %>%
    filter(TF %in% active_tfs) %>% # Only active TFs
    rename(TF = Motif, target = DORC, global_score = Score) %>%
    left_join(condition_tf_activity, by = "TF") %>%
    filter(!is.na(tf_avg_activity)) %>%
    mutate(
      celltype_refined = condition_name,
      tf_avg_activity = pmax(tf_avg_activity, 0),
      link_activity = tf_avg_activity * log2(global_score + 1)
    ) %>%
    # Filter for top links
    filter(link_activity > quantile(link_activity, quantile_thresh)) %>%
    arrange(desc(link_activity))
  
  return(condition_network)
}

# 6. Generate Networks ---------------------------------------------------------

message("Constructing stage-specific networks...")

# Generate networks for key stages
network_rs_early <- create_condition_network("Round_Spermatids_Early", figR.d_filtered, tf_activity_stats, TF_ACTIVITY_THRESHOLD, NETWORK_QUANTILE)
network_es_early <- create_condition_network("Elongating_Spermatids_Early", figR.d_filtered, tf_activity_stats, TF_ACTIVITY_THRESHOLD, NETWORK_QUANTILE)
network_es_late  <- create_condition_network("Elongating_Spermatids_Late", figR.d_filtered, tf_activity_stats, TF_ACTIVITY_THRESHOLD, NETWORK_QUANTILE)

# Combine into master dataframe
master_network <- bind_rows(network_rs_early, network_es_early, network_es_late)

# Save Master Network
saveRDS(master_network, file = glue("{OUTPUT_DIR}/stage_specific_networks.rds"))

# 7. Visualization Function ----------------------------------------------------

plot_network <- function(network_data, 
                         condition_name, 
                         top_n = 150, 
                         min_targets = 3, 
                         layout = "fr") {
  
  # Filter data for condition
  cond_data <- network_data %>%
    filter(celltype_refined == condition_name)
  
  # Filter TFs based on minimum connectivity (Out-Degree)
  valid_tfs <- cond_data %>%
    group_by(TF) %>%
    summarise(n_targets = n_distinct(target)) %>%
    filter(n_targets >= min_targets) %>%
    pull(TF)
  
  # Select Top N edges among valid TFs
  data_to_plot <- cond_data %>%
    filter(TF %in% valid_tfs) %>%
    arrange(desc(link_activity)) %>%
    head(top_n)
  
  if (nrow(data_to_plot) == 0) return(NULL)
  
  # Create Graph Object
  g <- graph_from_data_frame(data_to_plot, directed = TRUE)
  
  # Node Attributes
  all_tfs <- unique(data_to_plot$TF)
  all_targets <- unique(data_to_plot$target)
  V(g)$type <- ifelse(V(g)$name %in% all_tfs, 
                      ifelse(V(g)$name %in% all_targets, "Both", "TF"), "DORC")
  V(g)$degree <- degree(g, mode = "all")
  
  set.seed(42) # Reproducibility
  
  # Plot
  p <- ggraph(g, layout = layout) +
    geom_edge_fan(aes(alpha = link_activity, width = link_activity, color = link_activity),
                  arrow = arrow(length = unit(3, 'mm'), type = "closed"),
                  end_cap = circle(3, 'mm'), start_cap = circle(2, 'mm')) +
    geom_node_point(aes(color = type, size = degree), alpha = 0.9, stroke = 0.5) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = 20,
                   bg.color = "white", bg.r = 0.15, fontface = "bold") +
    scale_edge_width_continuous(range = c(0.3, 2.5), guide = "none") +
    scale_edge_alpha_continuous(range = c(0.3, 1.0), guide = "none") +
    scale_edge_color_gradientn(
      colors = c("#2166ac", "#67a9cf", "#d1e5f0", "#fddbc7", "#ef8a62", "#b2182b"),
      name = "Link Activity"
    ) +
    scale_color_manual(values = c("TF" = "#E31A1C", "DORC" = "#1F78B4", "Both" = "#33A02C")) +
    scale_size_continuous(range = c(3, 10)) +
    theme_graph(base_size = 12) +
    labs(title = gsub("_", " ", condition_name),
         subtitle = glue("Top {nrow(data_to_plot)} Links | Min Targets ≥ {min_targets}")) +
    theme(legend.position = "bottom")
  
  return(p)
}

# 8. Generate and Save Plots ---------------------------------------------------

message("Generating network visualizations...")

plots_list <- list(
  "RS_Early" = plot_network(master_network, "Round_Spermatids_Early", TOP_N_EDGES_PLOT, MIN_TARGETS_FOR_PLOT),
  "ES_Early" = plot_network(master_network, "Elongating_Spermatids_Early", TOP_N_EDGES_PLOT, MIN_TARGETS_FOR_PLOT),
  "ES_Late"  = plot_network(master_network, "Elongating_Spermatids_Late", TOP_N_EDGES_PLOT, MIN_TARGETS_FOR_PLOT)
)

# Save plots
for (name in names(plots_list)) {
  if (!is.null(plots_list[[name]])) {
    ggsave(filename = glue("{OUTPUT_DIR}/plots/{name}_network.pdf"), 
           plot = plots_list[[name]], 
           width = 12, height = 10, dpi = 300, device = cairo_pdf)
  }
}

# 9. Network Statistics & Comparison -------------------------------------------

message("Calculating network statistics...")

# Basic stats
network_stats <- master_network %>%
  group_by(celltype_refined) %>%
  summarise(
    n_edges = n(),
    n_tfs = n_distinct(TF),
    n_targets = n_distinct(target),
    avg_link_activity = mean(link_activity),
    .groups = 'drop'
  )

# Identify Top TFs per condition
top_tfs <- master_network %>%
  group_by(celltype_refined, TF) %>%
  summarise(total_activity = sum(link_activity), .groups = 'drop') %>%
  group_by(celltype_refined) %>%
  top_n(10, total_activity) %>%
  arrange(celltype_refined, desc(total_activity))

# Save Stats
write.csv(network_stats, glue("{OUTPUT_DIR}/network_global_statistics.csv"), row.names = FALSE)
write.csv(top_tfs, glue("{OUTPUT_DIR}/top_active_tfs.csv"), row.names = FALSE)

message("Analysis complete. Results saved to ", OUTPUT_DIR)
