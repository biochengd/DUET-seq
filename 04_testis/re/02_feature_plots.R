# ==============================================================================
# DUET-SEQ Analysis Pipeline: RE Population Feature Visualization
# ==============================================================================
# Description: 
# This script generates FeaturePlots for specific marker genes associated with 
# Round Spermatids (RS) and Elongating Spermatids (ES) to visualize their 
# expression patterns on the UMAP embedding.
#
# Steps:
# 1. Load the processed RE subset Seurat object.
# 2. Define marker genes of interest.
# 3. Generate customized FeaturePlots.
# ==============================================================================

# 1. Library Loading
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(qs)
  library(glue)
  library(patchwork)
})

# 2. Directory Setup & Data Loading
# ==============================================================================
# Set working directory relative to project root
# setwd("./src/rna/seurat_signac_pipeline2/04_re/01_preprocessing_qc_clustering")

data_dir <- "../data"
outs_dir <- "../outs"

if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(outs_dir)) dir.create(outs_dir, recursive = TRUE)

message("Loading RE Seurat object...")
# Load the processed object (ensure the file exists from previous steps)
seu_re <- qread(file = glue("{data_dir}/seu_re.qs"))

# 3. Visualization Configuration
# ==============================================================================
DefaultAssay(seu_re) <- "SCT"
reduction_used <- "umap.rna.dims_1_35"

# Define markers to plot
markers_filtered <- c("Cd63", "Gopc", "Sun5", "Spag6", "Tnp1", "Prm2", "Gapdhs", "Adam3")

# Define a clean theme for publication-quality plots
clean_theme <- theme(
  plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 10),
  legend.position = "right",
  aspect.ratio = 1,
  # Remove axis lines, text, and ticks
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  # Remove axis titles
  axis.title.x = element_blank(),
  axis.title.y = element_blank()
)

# 4. Generate Feature Plots
# ==============================================================================
message("Generating FeaturePlots...")

featureplot1 <- FeaturePlot(
  seu_re, 
  features = markers_filtered, 
  pt.size = 0.5, 
  order = TRUE, 
  cols = c("lightgrey", "#DE1F1F"), 
  raster = FALSE, 
  reduction = reduction_used
) & clean_theme

# 5. Save Results
# ==============================================================================
message("Saving plot...")

ggsave(
  filename = glue("{outs_dir}/featureplot_testis_re_filtered.pdf"), 
  plot = featureplot1, 
  width = 20, 
  height = 10 # Adjusted height to likely fit 2 rows of 4 plots better
)

message("Visualization completed successfully.")
