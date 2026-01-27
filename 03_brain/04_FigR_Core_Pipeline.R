#!/usr/bin/env Rscript

# ======== 1. CONFIGURATION ========
WORKING_DIR  <- "."                       # Project root
INPUT_FILE   <- "data/brain_multiome_cleaned.rds"
OUTPUT_DIR   <- "results/08_FigR_Analysis"
N_CORES      <- 2                         # Number of threads
GENOME       <- "mm10"                    # "mm10" or "hg38"

# ======== 2. SETUP ========
setwd(WORKING_DIR)
dir.create(file.path(OUTPUT_DIR, "Data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "QC"), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(FigR)
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(Matrix)
  library(parallel)
  library(doParallel)
  library(chromVAR)
  library(SummarizedExperiment)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(FNN)
})

registerDoParallel(cores = N_CORES)
options(future.globals.maxSize = 16000 * 1024^2)

cat("[1/6] Loading Data...\n")
if (!file.exists(INPUT_FILE)) stop("Input file not found.")
seu <- readRDS(INPUT_FILE)
cat(sprintf("   - Loaded %d cells.\n", ncol(seu)))

# ======== 3. PREPROCESSING ========
cat("[2/6] Preparing Multiome Data...\n")

# --- 3.1 Extract Counts ---
# Ensure common cells
common_cells <- intersect(colnames(seu[["RNA"]]), colnames(seu[["ATAC"]]))
rna_mat <- GetAssayData(seu, assay = "RNA", layer = "counts")[, common_cells]
atac_mat <- GetAssayData(seu, assay = "ATAC", layer = "counts")[, common_cells]

# Filter zero-expression genes
rna_mat <- rna_mat[Matrix::rowSums(rna_mat) > 0, ]

# --- 3.2 Chromosome Filtering ---
# Keep only standard chromosomes to avoid BSgenome errors
cat("   - Filtering non-standard chromosomes...\n")
peak_ranges <- StringToGRanges(rownames(atac_mat))
main_chrs <- c(paste0("chr", 1:19), "chrX", "chrY")
keep_peaks <- as.character(seqnames(peak_ranges)) %in% main_chrs
atac_mat_clean <- atac_mat[keep_peaks, ]

# --- 3.3 Create SummarizedExperiment ---
cat("   - Creating SummarizedExperiment object...\n")
atac_ranges <- StringToGRanges(rownames(atac_mat_clean))
seqlevelsStyle(atac_ranges) <- "UCSC"

ATAC.se <- SummarizedExperiment(
  assays = list(counts = atac_mat_clean),
  rowRanges = atac_ranges,
  colData = seu@meta.data[common_cells, ]
)

# Add UMAP coordinates if available (prioritize ATAC > RNA > WNN)
umap_red <- if("atac.umap" %in% names(seu@reductions)) "atac.umap" else "umap"
if(umap_red %in% names(seu@reductions)) {
  coords <- Embeddings(seu, reduction = umap_red)[common_cells, ]
  colData(ATAC.se)$UMAP1 <- coords[, 1]
  colData(ATAC.se)$UMAP2 <- coords[, 2]
}

# ======== 4. KNN COMPUTATION (LSI-BASED) ========
cat("[3/6] Computing KNN (using LSI)...\n")
# Use LSI for scalability (faster than cisTopic)
if ("lsi" %in% names(seu@reductions)) {
  lsi_coords <- Embeddings(seu, reduction = "lsi")[common_cells, 1:20]
  set.seed(123)
  cellkNN <- get.knn(lsi_coords, k = 30)$nn.index
  rownames(cellkNN) <- common_cells
} else {
  warning("LSI reduction not found. Smoothing might be skipped.")
  cellkNN <- NULL
}

# ======== 5. PEAK-GENE CORRELATION ========
cat("[4/6] Running Peak-Gene Correlation...\n")
cisCorr <- runGenePeakcorr(
  ATAC.se = ATAC.se,
  RNAmat = rna_mat,
  genome = GENOME,
  nCores = N_CORES,
  n_bg = 100,
  p.cut = NULL # Keep all for now, filter later
)

write.csv(cisCorr, file.path(OUTPUT_DIR, "Data/cisCorr_raw.csv"), row.names = FALSE)

# Filter for significant links (DORCs)
cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)
cat(sprintf("   - Found %d significant links.\n", nrow(cisCorr.filt)))

# Identify DORC genes (cutoff: >= 5 peaks)
dorcGenes <- dorcJPlot(cisCorr.filt, cutoff = 5, returnGeneList = TRUE)
cat(sprintf("   - Identified %d DORC genes.\n", length(dorcGenes)))

# ======== 6. SCORING & SMOOTHING ========
cat("[5/6] Calculating & Smoothing DORC Scores...\n")

# Calculate raw scores
dorcMat <- getDORCScores(ATAC.se, dorcTab = cisCorr.filt, geneList = dorcGenes, nCores = N_CORES)

# Smooth scores if KNN is available
if (!is.null(cellkNN)) {
  cat("   - Smoothing data using KNN...\n")
  dorcMat.s <- smoothScoresNN(NNmat = cellkNN, mat = dorcMat, nCores = N_CORES)
  rna_mat.s <- smoothScoresNN(NNmat = cellkNN, mat = rna_mat, nCores = N_CORES)
  
  saveRDS(dorcMat.s, file.path(OUTPUT_DIR, "Data/dorcMat_smoothed.rds"))
  saveRDS(rna_mat.s, file.path(OUTPUT_DIR, "Data/rnaMat_smoothed.rds"))
} else {
  dorcMat.s <- dorcMat
  rna_mat.s <- rna_mat
}

# ======== 7. TF-DORC NETWORK INFERENCE ========
cat("[6/6] Inferring TF-DORC Regulatory Network...\n")

figR.d <- runFigRGRN(
  ATAC.se = ATAC.se,
  dorcTab = cisCorr.filt,
  genome = GENOME,
  dorcMat = dorcMat.s,
  rnaMat = rna_mat.s,
  nCores = N_CORES
)

# Save Final Results
saveRDS(figR.d, file.path(OUTPUT_DIR, "Data/figR_results.rds"))
write.csv(figR.d, file.path(OUTPUT_DIR, "Data/figR_TF_DORC_table.csv"), row.names = FALSE)

cat("\n=== Core Analysis Complete ===\n")
