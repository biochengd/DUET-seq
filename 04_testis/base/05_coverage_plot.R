# ==============================================================================
# DUET-SEQ Analysis Pipeline: Genomic Coverage Plots
# ==============================================================================
# Description: 
# This script generates CoveragePlots (genome tracks) for specific marker genes 
# to visualize chromatin accessibility across different cell types.
#
# Steps:
# 1. Load the processed Seurat object.
# 2. Define cell type levels and color palette.
# 3. Generate and save coverage plots for a list of genes of interest.
# ==============================================================================

# 1. Library Loading
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(dplyr)
  library(glue)
  library(qs)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(EnsDb.Mmusculus.v79)
})

# 2. Directory Setup & Data Loading
# ==============================================================================
# Set working directory relative to project root
# setwd("./src/rna/seurat_signac_pipeline2/01_overall") 

data_dir <- "./data_overall"
outs_dir <- "./outs_overall"

if (!dir.exists(outs_dir)) dir.create(outs_dir, recursive = TRUE)

message("Loading Seurat object...")
# Load the object saved in the previous step (e.g., the one with links or the leiden clustered one)
seu_obj <- qread(file.path(data_dir, "seu_para1_remove2_33_9_15_remove7_19_leiden.qs"))

# 3. Setup Annotations and Colors
# ==============================================================================
DefaultAssay(seu_obj) <- 'ATAC'

# Define standard cell type levels
celltype_levels <- c(
  "Spermatogonia_SPG", "Spermatocytes_SPC", "Round_spermatids", 
  "Elongating_spermatids", "Sertoli_cells", "Leydig_cells", 
  "Stromal_cells", "Peritubular_myoid_cells", "Macrophage", "Endothelial_cells"
)

# Ensure metadata matches levels (handle renaming if necessary)
# Note: Ensure 'Telocyte' is renamed to 'Stromal_cells' if that is the intention
if ("Telocyte" %in% unique(seu_obj$celltype)) {
    seu_obj$celltype <- gsub("Telocyte", "Stromal_cells", seu_obj$celltype)
}

# Set factor levels
seu_obj$celltype <- factor(seu_obj$celltype, levels = celltype_levels)
Idents(seu_obj) <- "celltype"

# Define color palette
cell_type_colors <- c(
  "Spermatogonia_SPG" = '#FB8D3C', 
  "Spermatocytes_SPC" = '#d6604d',
  "Round_spermatids" = "#D082AF",
  "Elongating_spermatids" = "#7B599C", 
  "Sertoli_cells" = "#00ACB3",
  "Leydig_cells" = "#65AB53",
  "Stromal_cells" = "#4667A8",
  "Peritubular_myoid_cells" = "#9CC0DD", 
  "Macrophage" = "#85BFB8",
  "Endothelial_cells" = '#820610'
)

# 4. Generate Coverage Plots
# ==============================================================================
message("Generating coverage plots...")

genes_to_plot <- c('Uchl1', 'Sycp2', 'Piwil1', 'Acrv1', 'Prm2', 'Clu')

for (gene in genes_to_plot) {
  # Use tryCatch to ensure loop continues even if one gene fails (e.g., symbol not found)
  tryCatch({
    message(glue("Plotting {gene}..."))
    
    cov_plot <- CoveragePlot(
      object = seu_obj,
      region = gene,
      annotation = TRUE,
      features = gene,
      expression.assay = "SCT",
      extend.upstream = 10000,
      extend.downstream = 10000,
      assembly = "mm10",
      peaks = TRUE
    ) + scale_fill_manual(values = cell_type_colors)
    
    # Save plot
    ggsave(
      filename = glue("{outs_dir}/{gene}_combined_track_plot.pdf"), 
      plot = cov_plot, 
      width = 10, 
      height = 10
    )
    
    message(glue("Successfully saved plot for {gene}"))
    
  }, error = function(e) {
    message(glue("Error plotting {gene}: {e$message}"))
  })
}

message("Visualization pipeline completed.")
