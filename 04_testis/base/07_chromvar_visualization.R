# ==============================================================================
# DUET-SEQ Analysis Pipeline: ChromVAR Motif Activity Visualization
# ==============================================================================
# Description: 
# This script visualizes the activity of transcription factor (TF) motifs 
# inferred by ChromVAR alongside the expression of their corresponding TFs.
# It generates:
# 1. Gene Expression Density Plots (Nebulosa).
# 2. ChromVAR Motif Activity Density Plots.
# 3. Gene Activity Score Density Plots (if available).
# 4. Motif PWM Logo Plots.
#
# Steps:
# 1. Load the processed Seurat object (with ChromVAR assay) and motif mapping data.
# 2. Define plotting function for consistent styling.
# 3. Loop through specific TFs of interest and generate plots.
# ==============================================================================

# 1. Library Loading
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(dplyr)
  library(qs)
  library(glue)
  library(Nebulosa) # For plot_density
  library(TFBSTools)
  library(JASPAR2020)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

# 2. Directory Setup & Data Loading
# ==============================================================================
# Set working directory relative to project root
# setwd("./src/rna/seurat_signac_pipeline2/01_overall") 

data_dir <- "./data_overall"
outs_dir <- "./outs_overall/chromvar_plots"

if (!dir.exists(outs_dir)) dir.create(outs_dir, recursive = TRUE)

message("Loading data...")
# Load Seurat object (must have 'chromvar_sertoli' or relevant assay computed)
seu_obj <- qread(file = glue("{data_dir}/seu_sertoli_no_r3_sub_modified_atac.qs"))

# Load differential motif results or mapping table (motif_id <-> gene_symbol)
# Assuming 'da_chromvar_common' dataframe exists from previous analysis step
# which links motif IDs to gene symbols and filters for significant ones.
da_chromvar_common <- qread(file = glue("{data_dir}/da_chromvar_common.qs")) 

# Load motif metadata if needed for MotifPlot (usually stored in the object or a separate DF)
# If stored separately:
# motif_df <- qread(file = glue("{data_dir}/motif_df_sertoli.qs"))

# Ensure correct assays are set/available
# Note: 'chromvar_sertoli' is the assay name used in the original script
chromvar_assay_name <- "chromvar_sertoli" 

# 3. Visualization Configuration
# ==============================================================================

# Define a consistent theme for density plots
density_theme <- theme(
  legend.frame = element_rect(colour = "black"),
  legend.ticks = element_line(colour = "black", linewidth = 0),
  legend.key.width = unit(0.3, "cm"),
  legend.key.height = unit(0.8, "cm"),
  legend.title = element_text(color = 'black', face = "bold", size = 8),
  panel.border = element_rect(colour = "black", fill = NA),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  axis.title = element_blank()
)

# 4. Plotting Loop
# ==============================================================================
message("Generating ChromVAR visualization plots...")

# Extract unique TFs to plot from the mapping dataframe
tf_symbols <- unique(da_chromvar_common$gene_symbol)
# Filter out NAs or empty strings
tf_symbols <- tf_symbols[!is.na(tf_symbols) & tf_symbols != ""]

# Optional: Subset for testing
# tf_symbols <- head(tf_symbols, 5) 

for (gene_aimed in tf_symbols) {
  
  message(glue("Processing TF: {gene_aimed}"))
  
  # --- 1. Gene Expression (SCT) ---
  tryCatch({
    DefaultAssay(seu_obj) <- "SCT"
    
    if (gene_aimed %in% rownames(seu_obj)) {
      p_expr <- plot_density(
        seu_obj, 
        features = gene_aimed, 
        reduction = "umap.rna.dims_1_35", # Adjust reduction name as needed
        pal = 'viridis', 
        raster = TRUE, 
        size = 1
      ) + density_theme + ylab("Gene Expression")
      
      ggsave(
        filename = glue("{outs_dir}/density_expression_{gene_aimed}.pdf"), 
        plot = p_expr, 
        width = 5, height = 5
      )
    } else {
      message(glue("  -> Gene {gene_aimed} not found in SCT assay."))
    }
  }, error = function(e) { message(glue("  -> Error plotting expression: {e$message}")) })
  
  
  # --- 2. ChromVAR Motif Activity ---
  tryCatch({
    # Find the corresponding motif ID
    # Assuming 'da_chromvar_common' has columns 'gene_symbol' and 'motif_id' (or 'gene' as motif ID)
    # Adjust column names based on your actual object structure
    motif_id <- da_chromvar_common$motif_id[da_chromvar_common$gene_symbol == gene_aimed][1]
    
    if (!is.na(motif_id) && motif_id %in% rownames(seu_obj[[chromvar_assay_name]])) {
      DefaultAssay(seu_obj) <- chromvar_assay_name
      
      p_activity <- plot_density(
        seu_obj, 
        features = motif_id, 
        reduction = "umap.atac.dims_2_12", # Adjust reduction name
        pal = 'magma', 
        raster = TRUE, 
        size = 1
      ) + density_theme + ylab("Motif Activity")
      
      ggsave(
        filename = glue("{outs_dir}/density_chromvar_{gene_aimed}.pdf"), 
        plot = p_activity, 
        width = 5, height = 5
      )
      
      # --- 3. Motif Logo Plot ---
      # Requires the ATAC assay to have motif information added
      p_logo <- MotifPlot(
        object = seu_obj,
        motifs = motif_id,
        assay = 'ATAC'
      ) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
      ) +
      ggtitle(glue("{gene_aimed} ({motif_id})"))
      
      ggsave(
        filename = glue("{outs_dir}/logo_motif_{gene_aimed}_{motif_id}.pdf"), 
        plot = p_logo, 
        width = 5, height = 3
      )
      
    } else {
      message(glue("  -> Motif ID for {gene_aimed} not found or not in assay."))
    }
  }, error = function(e) { message(glue("  -> Error plotting motif: {e$message}")) })
  
  
  # --- 4. Gene Activity Score (Optional) ---
  # Only if 'ACTIVITY' assay exists and gene is present
  if ("ACTIVITY" %in% names(seu_obj@assays)) {
    tryCatch({
      DefaultAssay(seu_obj) <- "ACTIVITY"
      if (gene_aimed %in% rownames(seu_obj)) {
        p_gene_act <- plot_density(
          seu_obj, 
          features = gene_aimed,
          reduction = "umap.atac.dims_2_12",
          pal = 'inferno', 
          raster = TRUE, 
          size = 1
        ) + density_theme + ylab("Gene Activity Score")
        
        ggsave(
          filename = glue("{outs_dir}/density_gene_activity_{gene_aimed}.pdf"), 
          plot = p_gene_act, 
          width = 5, height = 5
        )
      }
    }, error = function(e) { message(glue("  -> Error plotting gene activity: {e$message}")) })
  }
}

message("ChromVAR visualization pipeline completed.")
