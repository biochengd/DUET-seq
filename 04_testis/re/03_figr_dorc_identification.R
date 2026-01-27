# ==============================================================================
# DUET-SEQ Analysis Pipeline: FigR - DORC Identification (Step 1)
# ==============================================================================
# Description: 
# This script is the first part of the FigR workflow. It identifies Domains of 
# Regulatory Chromatin (DORCs) by linking ATAC-seq peaks to gene expression.
#
# Steps:
# 1. Environment setup and parallel processing configuration.
# 2. Load processed Seurat object.
# 3. Prepare input data for FigR (SummarizedExperiment, RNA matrix, KNN).
# 4. Perform Peak-Gene correlation analysis.
# 5. Identify DORCs and save intermediate results for the GRN step.
#
# Next Step: 
# After running this script, check the "DORC_JPlot.pdf" output. If the number 
# of DORCs is satisfactory, proceed to '12_figr_grn_analysis.R'.
# ==============================================================================

# 1. Library Loading & Environment
# ==============================================================================
suppressPackageStartupMessages({
  library(FigR)
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(GenomicRanges)
  library(Matrix)
  library(parallel)
  library(doParallel)
  library(SummarizedExperiment)
  library(FNN)
  library(qs)
  library(glue)
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

# Directory Setup (Relative Paths)
# setwd("./src/rna/seurat_signac_pipeline2/04_re/04_figr_overall") 

data_dir <- "../data" 
figr_out_dir <- "./08_FigR_Analysis"
dirs_to_create <- c(
  figr_out_dir,
  file.path(figr_out_dir, "Peak_Gene_Links"),
  file.path(figr_out_dir, "DORC_Analysis"), 
  file.path(figr_out_dir, "Data")
)

for(dir in dirs_to_create) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# Parallel Processing Setup
registerDoParallel(cores = 2)
options(future.globals.maxSize = 16000 * 1024^4)

# 2. Load Data
# ==============================================================================
message("Loading Seurat object...")
# Assuming 'seu_re.qs' (or specific subset) is the input
seu <- qread(file = glue("{data_dir}/seu_re.qs")) # Updated filename reference

# Ensure correct cell type column is used
if ("celltype_re_lvl2" %in% colnames(seu@meta.data)) {
    seu$celltype <- seu$celltype_re_lvl2
}

message(paste0("Data loaded: ", ncol(seu), " cells"))

# Ensure annotations are present
if(is.null(Annotation(seu[["ATAC"]]))) {
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations) <- "UCSC"
  Annotation(seu[["ATAC"]]) <- annotations
}

# 3. Prepare FigR Inputs
# ==============================================================================
message("Preparing inputs for FigR...")

# Extract Counts
DefaultAssay(seu) <- "SCT"
rna_counts <- GetAssayData(seu, assay = "SCT", layer = "counts")

DefaultAssay(seu) <- "ATAC"
atac_counts <- GetAssayData(seu, assay = "ATAC", layer = "counts")

# Align Cells
common_cells <- intersect(colnames(rna_counts), colnames(atac_counts))
rna_counts <- rna_counts[, common_cells]
atac_counts <- atac_counts[, common_cells]

# Filter non-expressed genes
rna_counts <- rna_counts[Matrix::rowSums(rna_counts) != 0, ]

# Chromosome Filtering (Standard chromosomes only)
message("Filtering peaks for standard chromosomes...")
main_chromosomes <- c(paste0("chr", 1:19), "chrX", "chrY", "chrM")
peak_ranges <- StringToGRanges(rownames(atac_counts))
main_chr_mask <- as.character(seqnames(peak_ranges)) %in% main_chromosomes
atac_counts_clean <- atac_counts[main_chr_mask, ]

message(paste0("Peaks remaining: ", nrow(atac_counts_clean)))

# Create SummarizedExperiment for ATAC
atac_ranges_clean <- StringToGRanges(rownames(atac_counts_clean))
seqlevelsStyle(atac_ranges_clean) <- "UCSC"
ATAC.se <- SummarizedExperiment(
  assays = list(counts = atac_counts_clean),
  rowRanges = atac_ranges_clean,
  colData = seu@meta.data[common_cells, ]
)

RNAmat <- rna_counts

# Calculate Cell kNN (using LSI)
message("Calculating cell kNN matrix...")
if ("lsi" %in% names(seu@reductions)) {
  lsi_coords <- Embeddings(seu, reduction = "lsi")[common_cells, ]
  set.seed(123)
  # Use top 20 LSI components
  cellkNN <- get.knn(lsi_coords[, 1:min(20, ncol(lsi_coords))], k = 30)$nn.index
  rownames(cellkNN) <- common_cells
} else {
  stop("LSI reduction not found in Seurat object.")
}

# 4. Peak-Gene Correlation Analysis
# ==============================================================================
message("Running Peak-Gene correlation analysis...")

cisCorr_path <- file.path(figr_out_dir, "Peak_Gene_Links/figr_peak_gene_links.csv")

if (!file.exists(cisCorr_path)) {
  cisCorr <- FigR::runGenePeakcorr(
    ATAC.se = ATAC.se,
    RNAmat = RNAmat,
    genome = "mm10", # Ensure genome matches your data
    nCores = 2,
    p.cut = NULL,
    n_bg = 100
  )
  write.csv(cisCorr, cisCorr_path, row.names = FALSE)
} else {
  message("Loading existing Peak-Gene links...")
  cisCorr <- read.csv(cisCorr_path)
}

# 5. DORC Identification
# ==============================================================================
message("Identifying DORCs...")

# Filter significant links (Z-test p-value <= 0.05)
cisCorr.filt <- cisCorr %>% dplyr::filter(pvalZ <= 0.05)
message(paste0("Significant Peak-Gene links: ", nrow(cisCorr.filt)))

# Generate DORC Plot and Extract Genes
pdf(file.path(figr_out_dir, "DORC_Analysis/DORC_JPlot.pdf"), width = 10, height = 8)
dorcGenes <- dorcJPlot(
  dorcTab = cisCorr.filt,
  cutoff = 4,        # Minimum number of significant peaks per gene to call a DORC
  labelTop = 20,
  returnGeneList = TRUE,
  force = 2
)
dev.off()

message(paste0("Identified DORC genes: ", length(dorcGenes)))

# 6. Save Intermediate Data for Step 2
# ==============================================================================
message("Saving intermediate data for GRN analysis...")

save_path <- file.path(figr_out_dir, "Data")
saveRDS(ATAC.se, file.path(save_path, "ATAC_se.rds"))
saveRDS(RNAmat, file.path(save_path, "RNAmat.rds"))
saveRDS(cellkNN, file.path(save_path, "cellkNN.rds"))
saveRDS(cisCorr.filt, file.path(save_path, "cisCorr_filt.rds"))
saveRDS(dorcGenes, file.path(save_path, "dorcGenes.rds"))

message("Step 1 complete. Proceed to '12_figr_grn_analysis.R'.")
