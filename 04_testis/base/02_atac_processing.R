# ==============================================================================
# DUET-SEQ Analysis Pipeline: ATAC Processing & Clustering
# ==============================================================================
# Description: 
# This script performs normalization, dimensionality reduction (LSI), Harmony 
# integration, and clustering specifically for the ATAC modality on the 
# filtered Multiome data.
#
# Steps:
# 1. Load the filtered Seurat object list (from RNA QC step).
# 2. ATAC Normalization (TF-IDF) and Feature Selection.
# 3. Dimensionality Reduction (SVD/LSI).
# 4. Batch Correction (Harmony) on ATAC reductions.
# 5. Clustering and UMAP embedding generation across multiple resolutions/dims.
# ==============================================================================

# 1. Library Loading
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(harmony)
  library(qs)
  library(glue)
  library(ggplot2)
  # Genome annotations
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

# 2. Directory Setup & Data Loading
# ==============================================================================
# Set working directory relative to project root
# setwd("./src/rna/seurat_signac_pipeline2") 

# Define input/output directories
data_dir <- "./data"
outs_dir <- "./outs"

if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(outs_dir)) dir.create(outs_dir, recursive = TRUE)

message("Loading filtered Seurat object list...")
# Assuming 'seu_filtered_list.qs' contains the list of Seurat objects after RNA QC
# If continuing from the single object workflow, load that object instead.
# Here we follow the logic of processing a list (e.g., parameter sets).
seu_filtered_list <- qread(file = glue("{data_dir}/seu_filtered_list.qs"))

# 3. Parallel Processing Setup
# ==============================================================================
library(future)
plan("multisession", workers = 2)
options(future.globals.maxSize = 80000 * 1024^4)

# 4. ATAC Processing Pipeline
# ==============================================================================
message("Starting ATAC processing pipeline...")

para_names <- names(seu_filtered_list)
sig_filtered_list2 <- list()

for (i in seq_along(seu_filtered_list)) {
    current_seu <- seu_filtered_list[[i]]
    dataset_name <- para_names[i]
    
    message(paste0("Processing dataset: ", dataset_name, " (", i, "/", length(seu_filtered_list), ")"))
    
    DefaultAssay(current_seu) <- "ATAC"
    group_by_vars_col <- c("sample_names", "group")

    # 4.1 Normalization and LSI
    # --------------------------------------------------------------------------
    message(" -> Running TF-IDF, Feature Selection, and SVD...")
    current_seu <- current_seu |>
        RunTFIDF(method = 1) |>
        FindTopFeatures(min.cutoff = 'q5') |>
        RunSVD(n = 50)

    # 4.2 Harmony Integration
    # --------------------------------------------------------------------------
    message(" -> Running Harmony integration on ATAC LSI...")
    current_seu <- current_seu |>
    RunHarmony(
            group.by.vars = group_by_vars_col,
            reduction.use = "lsi",
            reduction.save = "harmony.atac",
            assay.use = "ATAC", 
            dims.use = 2:50, # Skip 1st LSI component as it often correlates with sequencing depth
            project.dim = FALSE,
            verbose = FALSE
          )

    # 4.3 Clustering and UMAP
    # --------------------------------------------------------------------------
    # Define range of dimensions to test
    atac_dims_options <- list(
        "dims_2_8" = 2:8, "dims_2_10" = 2:10, "dims_2_12" = 2:12,
        "dims_2_14" = 2:14, "dims_2_16" = 2:16, "dims_2_18" = 2:18,
        "dims_2_20" = 2:20, "dims_2_25" = 2:25, "dims_2_30" = 2:30
    )
    
    res_opts <- c(0.4, 0.6, 0.8, 1, 2)
    
    message(" -> Running Neighbors, Clusters, and UMAP for multiple dimension sets...")
    
    for(dim_name in names(atac_dims_options)) {
      dims_to_use <- atac_dims_options[[dim_name]]

      # Define graph names for this dimension set
      graph_nn_name <- paste0("atac_nn_", dim_name)
      graph_snn_name <- paste0("atac_snn_", dim_name)

      current_seu <- current_seu |>
          FindNeighbors(
              reduction = "harmony.atac",
              dims = dims_to_use,
              graph.name = c(graph_nn_name, graph_snn_name),
              verbose = FALSE
          ) |>
          FindClusters(
              graph.name = graph_snn_name,
              resolution = res_opts,
              algorithm = 3, # SLM algorithm
              verbose = FALSE
          ) |>
          RunUMAP(
            reduction = "harmony.atac", 
            dims = dims_to_use, 
            reduction.name = paste0("umap.atac.", dim_name), 
            reduction.key = paste0("atacUMAP", gsub("dims_", "", dim_name), "_"),
            n.neighbors = 50,
            min.dist = 0.12,
            spread = 1,
            verbose = FALSE
          )
    }
    
    sig_filtered_list2[[dataset_name]] <- current_seu
    
    # Clean up memory
    gc()
}

# 5. Save Results
# ==============================================================================
message("Saving processed ATAC objects...")
qsave(sig_filtered_list2, file = glue("{data_dir}/sig_filtered_list2.qs"))

message("Pipeline completed successfully.")
