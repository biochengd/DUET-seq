
# ==============================================================================
# DUET-SEQ Analysis Pipeline: Preprocessing, QC, and Clustering
# ==============================================================================
# Description: 
# This script performs quality control (QC), filtering, and doublet removal 
# for single-cell Multiome (RNA + ATAC) data.
#
# Steps:
# 1. Load merged Seurat object.
# 2. Metadata formatting and grouping.
# 3. RNA modality QC and filtering.
# 4. ATAC modality QC and filtering (using specific parameter set).
# 5. Doublet detection and removal using scDblFinder.
# 6. Normalization, Dimensionality Reduction, and Clustering on clean data.
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
  library(scDblFinder)
  library(dplyr)
  library(SingleCellExperiment)
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

message("Loading raw merged Seurat object...")
seu <- qread(file = glue("{data_dir}/seu_merged_v3.qs"))

# 3. Metadata Processing
# ==============================================================================
message("Processing metadata...")

DefaultAssay(seu) <- "RNA"

# Create group labels based on sample names
seu$group <- case_when(
  seu$sample_name %in% c("250704DDA_P7_1", "250704DDA_P7_2", "250704DDA_P7_3") ~ "P7",
  seu$sample_name %in% c("250605DDA_P14_1", "250605DDA_P14_2", "250605DDA_P14_3") ~ "P14",
  seu$sample_name %in% c("250605DDA_P21_1", "250605DDA_P21_2", "250605DDA_P21_3") ~ "P21",
  seu$sample_name %in% c("250606DDA_P28_1", "250606DDA_P28_2", "250606DDA_P28_3") ~ "P28",
  seu$sample_name %in% c("250709DDA_P35_1", "250709DDA_P35_2", "250709DDA_P35_3") ~ "P35",
  seu$sample_name %in% c("250606DDA_10W_1", "250606DDA_10W_2", "250606DDA_10W_3") ~ "10W",
  TRUE ~ "Unknown"
)

# Set factor levels for ordered visualization
my_levels <- c(
  "250704DDA_P7_1", "250704DDA_P7_2", "250704DDA_P7_3", 
  "250605DDA_P14_1", "250605DDA_P14_2", "250605DDA_P14_3", 
  "250605DDA_P21_1", "250605DDA_P21_2", "250605DDA_P21_3", 
  "250606DDA_P28_1", "250606DDA_P28_2", "250606DDA_P28_3", 
  "250709DDA_P35_1", "250709DDA_P35_2", "250709DDA_P35_3", 
  "250606DDA_10W_1", "250606DDA_10W_2", "250606DDA_10W_3"
)
seu$orig.ident <- factor(seu$orig.ident, levels = my_levels)

# Clean up unnecessary metadata columns
cols_to_remove <- c("duplicate", "chimeric", "unmapped", "lowmapq", "nonprimary", 
                    "is__cell_barcode", "excluded_reason", "DNase_sensitive_region_fragments", 
                    "on_target_fragments", "peak_region_cutsites", "nucleosome_percentile",
                    "mitochondrial", "enhancer_region_fragments", "promoter_region_fragments", 
                    "blacklist_region_fragments", "TSS.percentile", "peak_region_fragments")

seu@meta.data <- seu@meta.data[, !colnames(seu@meta.data) %in% cols_to_remove]

# Standardize sample name column
seu$sample_names <- seu$sample_name
seu$sample_name <- NULL

# 4. RNA QC Calculation & Basic Filtering
# ==============================================================================
message("Calculating RNA QC metrics...")

# Calculate percentages for mitochondria, ribosomes, and hemoglobin
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
seu[["percent.rb"]] <- PercentageFeatureSet(seu, pattern = "^Rps|^Rpl")
seu[["percent.hb"]] <- PercentageFeatureSet(seu, pattern = "^Hb[^p]")

# RNA Filtering
# Criteria: 
# - UMI Counts < 50,000
# - 300 < Genes < 9,000
# - Mitochondrial content < 25%
# - Ribosomal content < 10%
message("Applying RNA filtering thresholds...")
seu <- subset(seu, subset = nCount_RNA < 50000 & 
                            nFeature_RNA > 300 & nFeature_RNA < 9000 & 
                            percent.mt < 25 & percent.rb < 10)

# 5. ATAC QC Filtering
# ==============================================================================
message("Applying ATAC filtering thresholds (Parameter Set 1)...")

# Criteria (Para1):
# - Features > 400
# - UMI Counts < 100,000
# - TSS Enrichment > 2
# - Blacklist Ratio < 0.075
# - Nucleosome Signal < 1.5
# - FRiP > 0.15
seu <- subset(
  seu,
  subset = nFeature_ATAC > 400 &
           nCount_ATAC < 100000 &
           TSS.enrichment > 2 &
           blacklist_ratio < 0.075 &
           nucleosome_signal < 1.5 &
           FRiP > 0.15
)

message(paste0("Cells remaining after QC: ", ncol(seu)))

# 6. Doublet Detection (scDblFinder)
# ==============================================================================
message("Running Doublet Detection using scDblFinder...")

DefaultAssay(seu) <- "RNA"
# Ensure layers are joined before splitting
seu <- JoinLayers(seu)

# Split by sample to calculate doublets individually
seurat_list <- SplitObject(seu, split.by = "orig.ident")

seurat_list <- lapply(seurat_list, function(obj) {
  # Convert to SingleCellExperiment for scDblFinder
  sce <- as.SingleCellExperiment(obj, assay = "RNA")
  sce <- scDblFinder(sce)
  # Transfer scores back to Seurat object
  obj$scDblFinder.score <- sce$scDblFinder.score
  obj$scDblFinder.class <- sce$scDblFinder.class
  return(obj)
})

# Merge back into one object
seu_scored <- merge(x = seurat_list[[1]], y = seurat_list[-1])
DefaultAssay(seu_scored) <- "RNA"
seu_scored <- JoinLayers(seu_scored)

# Filter out doublets
message(paste0("Total doublets detected: ", sum(seu_scored$scDblFinder.class == "doublet")))
seu_singlets <- subset(seu_scored, subset = scDblFinder.class == "singlet")

message(paste0("Cells remaining after doublet removal: ", ncol(seu_singlets)))

# Save the clean singlet object
qsave(seu_singlets, file = glue("{data_dir}/seu_qc_singlets.qs"))

# 7. Normalization, Dimensionality Reduction & Clustering
# ==============================================================================
message("Running SCTransform, PCA, Harmony, and UMAP...")

# Setup parallel processing
plan("multisession", workers = 2)
options(future.globals.maxSize = 80000 * 1024^4)

DefaultAssay(seu_singlets) <- "RNA"
hvg_num <- 3000
group_by_vars_col <- c("sample_names", "group")

# Processing workflow
seu_processed <- seu_singlets |>
  SCTransform(vars.to.regress = c("percent.mt"),
              variable.features.n = hvg_num,
              verbose = FALSE) |>
  RunPCA(npcs = 50, reduction.name = "pca") |>
  RunHarmony(reduction = "pca", 
             assay.use = "SCT", 
             group.by.vars = group_by_vars_col, 
             reduction.save = "harmony.rna",
             verbose = FALSE)

# Clustering and UMAP
# Using dims 1:30 as the standard selection based on previous analysis
selected_dims <- 1:30
resolution_val <- 0.8

seu_processed <- seu_processed |>
  FindNeighbors(reduction = "harmony.rna", dims = selected_dims) |>
  FindClusters(resolution = resolution_val, algorithm = 3, verbose = FALSE) |>
  RunUMAP(reduction = "harmony.rna",
          dims = selected_dims,
          reduction.name = "umap.rna",
          n.neighbors = 50,
          min.dist = 0.3,
          spread = 1.5,
          verbose = FALSE)

# 8. Save Final Object
# ==============================================================================
message("Saving final processed object...")
qsave(seu_processed, file = glue("{data_dir}/seu_processed_final.qs"))

message("Pipeline completed successfully.")
