# ==============================================================================
# DUET-SEQ Analysis Pipeline: Spermatogonia (SPG) Sub-analysis
# ==============================================================================
# Description: 
# This script focuses on the Spermatogonia (SPG) population extracted from the 
# overall dataset. It performs:
# 1. Stage definition (Early vs Late).
# 2. Identification of marker genes for different timepoints and stages.
# 3. Visualization of SPG subtypes and developmental trajectories.
#
# Steps:
# 1. Load the SPG subset Seurat object.
# 2. Define developmental stages.
# 3. Find top marker genes for groups and stages.
# 4. Visualize embeddings (UMAP) and timepoint contributions.
# ==============================================================================

# 1. Library Loading
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(qs)
  library(glue)
  library(ggrastr)
})

# 2. Directory Setup & Data Loading
# ==============================================================================
# Set working directory relative to project root
# setwd("./src/rna/seurat_signac_pipeline2/02_spg") 

data_dir <- "./data_spg"
outs_dir <- "./outs_spg"

if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(outs_dir)) dir.create(outs_dir, recursive = TRUE)

message("Loading SPG Seurat object...")
# Assuming 'seu_spg.qs' is the subsetted object containing only Spermatogonia
seu_spg <- qread(file = glue("{data_dir}/seu_spg.qs"))

# Define Stage based on Group (Timepoint)
# P7, P14 -> Early; Others -> Late
seu_spg$stage <- ifelse(seu_spg$group %in% c("P7", "P14"), "Early", "Late")

# Save updated metadata
qsave(seu_spg, file = glue("{data_dir}/seu_spg_processed.qs"))

# 3. Marker Gene Identification
# ==============================================================================
message("Identifying marker genes...")

DefaultAssay(seu_spg) <- "SCT"

# Helper function to find and save markers
find_and_save_markers <- function(seu_obj, ident_col, file_prefix) {
  Idents(seu_obj) <- ident_col
  message(glue("Finding markers for: {ident_col}"))
  
  markers <- FindAllMarkers(
    object = seu_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = FALSE
  )
  
  # Select Top 20 by LogFC
  top20 <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log2FC)
  
  write.csv(top20, file = glue("{outs_dir}/top20_genes_{file_prefix}.csv"), row.names = FALSE)
  write.csv(markers, file = glue("{outs_dir}/all_markers_{file_prefix}.csv"), row.names = FALSE)
  message(glue("Saved markers for {file_prefix}"))
}

# 3.1 Markers per Timepoint (Group)
find_and_save_markers(seu_spg, "group", "per_timepoint")

# 3.2 Markers per Stage (Early/Late)
find_and_save_markers(seu_spg, "stage", "per_stage")

# 4. Visualization
# ==============================================================================
message("Generating visualizations...")

# Define SPG Subtypes Colors (Adjust if your subtypes differ)
spg_colors <- c(
  "SPG_SSC" = '#FB8D3C', 
  "SPG_Undifferentiated" = '#d6604d',
  "SPG_Differentiating" = "#D082AF",
  "SPG_Meiotic" = "#7B599C"
)

# 4.1 RNA UMAP by Cell Type
# Ensure reduction exists, otherwise RunUMAP first
if (!"umap.rna.dims_1_25" %in% names(seu_spg@reductions)) {
    message("Calculating RNA UMAP...")
    seu_spg <- RunUMAP(
        seu_spg, 
        reduction = "harmony.rna", 
        dims = 1:25, 
        reduction.name = "umap.rna.dims_1_25", 
        reduction.key = "rnaUMAP_",
        verbose = FALSE
    )
}

p_rna <- DimPlot(
  seu_spg, 
  reduction = "umap.rna.dims_1_25", 
  group.by = "celltype_spg_lvl2", # Assuming this column exists for subtypes
  cols = spg_colors,
  pt.size = 0.5,
  raster = FALSE,
  shuffle = TRUE
) +
  theme_classic() +
  ggtitle("SPG Subtypes (RNA)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(filename = glue("{outs_dir}/spg_ subtypes_rna.pdf"), plot = p_rna, width = 8, height = 6)

# 4.2 Timepoint Highlight Plots
message("Generating timepoint highlight plots...")

timepoints <- c("P7", "P14", "P21", "P28", "P35", "10W")
scanpy_palette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")

plot_list <- list()

for (i in seq_along(timepoints)) {
  tp <- timepoints[i]
  color <- scanpy_palette[i]
  
  # Highlight cells from current timepoint
  cells_to_highlight <- WhichCells(seu_spg, expression = group == tp)
  
  p <- DimPlot(
    seu_spg, 
    reduction = "umap.rna.dims_1_25", 
    cells.highlight = cells_to_highlight, 
    cols.highlight = color, 
    cols = "grey90", 
    pt.size = 0.1,
    sizes.highlight = 0.1,
    raster = TRUE
  ) +
    ggtitle(tp) +
    NoLegend() +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_list[[tp]] <- p
}

# Combine plots
combined_highlight <- wrap_plots(plot_list, ncol = 3)
ggsave(filename = glue("{outs_dir}/spg_timepoint_highlight.pdf"), plot = combined_highlight, width = 12, height = 8)

message("SPG analysis visualization completed.")
