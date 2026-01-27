#!/usr/bin/env Rscript

# ==============================================================================
# 1. Configuration & Environment
# ==============================================================================

# --- User Configuration ---
# Set your working directory
WORKING_DIR <- "." 
# Input Seurat object path (Cleaned data from previous steps)
INPUT_FILE  <- "03_Cleaned_Analysis/Results/brain_multiome_cleaned.rds"
# Output directory base
OUTPUT_BASE <- "06a_ChromVAR_Analysis"
# Number of threads for parallel processing
N_CORES     <- 2
# Genome assembly (Ensure matching BSgenome package is installed)
GENOME_PKG  <- "BSgenome.Mmusculus.UCSC.mm10"

# --- Setup ---
setwd(WORKING_DIR)

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(future)
  library(TFBSTools)
  library(JASPAR2020)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(motifmatchr)
  library(viridis)
  library(chromVAR)
  library(BiocParallel)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

# Parallelization settings
plan("multisession", workers = N_CORES)
options(future.globals.maxSize = 80000 * 1024^2)

# Create directory structure
sub_dirs <- c(
  "Correlation_Analysis", 
  "DEG_Analysis", 
  "TF_Motif_Analysis", 
  "Data", 
  "Intermediate"
)
dirs_to_create <- file.path(OUTPUT_BASE, sub_dirs)
sapply(c(OUTPUT_BASE, dirs_to_create), dir.create, showWarnings = FALSE, recursive = TRUE)

cat("=== scMultiome ChromVAR Analysis Pipeline Started ===\n")

# ==============================================================================
# 2. Color Palette Definitions
# ==============================================================================

# Palette 1: Correlation (White to Red)
# Used for: 0.75 - 1.0 correlation values
get_correlation_colors <- function() {
  colorRamp2(
    c(0.75, 0.85, 0.95, 1.0), 
    c("#FFFFFF", "#FDDBC7", "#F4A582", "#B2182B")
  )
}

# Palette 2: Expression (Divergent Multicolor)
# Used for: Scaled gene expression (-2 to 4)
get_expression_colors <- function() {
  colorRamp2(
    c(-2, -1, 0, 1, 2, 3, 4),
    c("#053061", "#2166AC", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#B2182B")
  )
}

# Palette 3: TF Activity (Professional Divergent)
# Used for: Scaled chromVAR deviation scores (-3 to 3)
get_tf_colors <- function() {
  colorRamp2(
    c(-3, -1.5, 0, 1.5, 3),
    c("#2166AC", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#B2182B")
  )
}

# Reference Cell Type Palettes (Seurat & Kelly)
# These are used to generate the standardized color dictionary
SEURAT_COLORS <- c('#c77fb0', '#59aa4e', '#914187', '#82ac33', '#8a619f',
                   '#3bb0be', '#4d599d', '#e17a21', '#7386be', '#d7a437',
                   '#704a90', '#417f34', '#c13f86', '#5eb17d', '#dd2d63',
                   '#50b4a9', '#d5333a', '#445298', '#cd802e', '#bcb335',
                   '#89802b', '#9e3f84', '#9eab5e')

KELLY_COLORS <- c('#af2337', '#ecc342', '#2967a0', '#2f3c28', '#96b437',
                  '#da93ab', '#e58932', '#80598f', '#7e331f', '#3b855a',
                  '#c0b286', '#a9c9ed', '#ec977f', '#848482', '#604628',
                  '#d26034', '#a64c6b', '#dbd245', '#eba83b', '#5d5092',
                  '#222222', '#f2f3f4')

# ==============================================================================
# 3. Data Loading & Preprocessing
# ==============================================================================

cat("[Step 1] Loading data...\n")

if (!file.exists(INPUT_FILE)) {
  stop(sprintf("Error: Input file '%s' not found.", INPUT_FILE))
}

seu <- readRDS(INPUT_FILE)
cat(sprintf("   -> Loaded Seurat object with %d cells.\n", ncol(seu)))

# --- Assign Cell Type Colors ---
cat("[Step 2] Generating standardized color dictionary...\n")

actual_celltypes <- sort(unique(na.omit(seu$celltype)))
n_celltypes <- length(actual_celltypes)

# Combine palettes to ensure enough colors
full_palette <- c(SEURAT_COLORS, KELLY_COLORS)
if (n_celltypes > length(full_palette)) {
  warning("More cell types than available colors. Recycling colors.")
  full_palette <- rep(full_palette, ceiling(n_celltypes / length(full_palette)))
}

# Map alphabetically to ensure consistency
color_dict <- full_palette[1:n_celltypes]
names(color_dict) <- actual_celltypes

# Define logical biological order for plots (Modify this list based on your biology)
target_order <- c(
  "Ex_L2/3_IT", "Ex_L4_IT", "Ex_L5_IT", "Ex_L5_PT", "Ex_L5/6_IT",
  "Ex_L6_IT", "Ex_L6_CT", "Ex_L6b", "Ex_PIR_Ndst4", "Ex_Hippocampus_DG", 
  "Ex_Cerebellum_Gabra6", "In_Pvalb", "In_Sst", "In_Vip", "In_Gabrg1/Glra3", 
  "In_Striatum_MSN", "OPCs", "Oligodendrocytes", "Astrocytes", "Microglia",
  "Endothelial", "VLMCs", "Choroid_Plexus"
)

# Intersect with actual data to get the final ordered list
ordered_celltypes <- target_order[target_order %in% actual_celltypes]
missing_types <- setdiff(actual_celltypes, ordered_celltypes)
if(length(missing_types) > 0) {
  ordered_celltypes <- c(ordered_celltypes, missing_types)
}

cat(sprintf("   -> Color dictionary created for %d cell types.\n", n_celltypes))

# ==============================================================================
# 4. ChromVAR Analysis (With Checkpointing)
# ==============================================================================

cat("[Step 3] Running/Loading chromVAR analysis...\n")

chromvar_rds_path <- file.path(OUTPUT_BASE, "Intermediate", "seu_with_chromvar.rds")

if(file.exists(chromvar_rds_path)) {
  cat("   -> Found cached chromVAR results. Loading...\n")
  seu <- readRDS(chromvar_rds_path)
  cat("   -> Loaded successfully.\n")
} else {
  cat("   -> No cache found. Starting calculation...\n")
  
  DefaultAssay(seu) <- "ATAC"
  
  # 1. Genome Compatibility Check
  cat("   -> Checking genome compatibility...\n")
  peaks_gr <- granges(seu)
  current_seqlevels <- unique(as.character(seqnames(peaks_gr)))
  bsgenome_seqlevels <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
  
  # Filter non-standard chromosomes/scaffolds if necessary
  valid_chromosomes <- intersect(c(paste0("chr", 1:19), "chrX", "chrY"), bsgenome_seqlevels)
  keep_peaks <- as.logical(seqnames(peaks_gr) %in% valid_chromosomes)
  
  if (sum(keep_peaks) < length(peaks_gr)) {
    cat(sprintf("   -> Filtering scaffolds. Keeping %d / %d peaks.\n", sum(keep_peaks), length(peaks_gr)))
    
    # Rebuild ChromatinAssay with filtered peaks
    new_atac <- CreateChromatinAssay(
      counts = GetAssayData(seu, "ATAC", "counts")[keep_peaks, ],
      ranges = peaks_gr[keep_peaks],
      data = GetAssayData(seu, "ATAC", "data")[keep_peaks, ]
    )
    seu[["ATAC"]] <- new_atac
  }
  
  # 2. Motif Analysis
  # Parallel config for chromVAR
  register(SerialParam()) 
  options(mc.cores = 1) # chromVAR is often more stable with 1 core inside Rstudio/Seurat wrappers
  
  tryCatch({
    # Get JASPAR motifs
    pfm <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
    
    # Add motifs
    cat("   -> Adding motifs to object...\n")
    seu <- AddMotifs(object = seu, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)
    
    # Run chromVAR
    cat("   -> Computing deviations (RunChromVAR)...\n")
    seu <- RunChromVAR(object = seu, genome = BSgenome.Mmusculus.UCSC.mm10)
    
    # Ensure RNA is normalized before saving
    DefaultAssay(seu) <- "RNA"
    if(!"data" %in% Layers(seu[["RNA"]])) {
      seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)
    }
    
    # Save checkpoint
    saveRDS(seu, chromvar_rds_path)
    cat(sprintf("   -> Saved chromVAR results to %s\n", chromvar_rds_path))
    
  }, error = function(e) {
    stop(paste("ChromVAR analysis failed:", e$message))
  })
}

# ==============================================================================
# 5. Correlation Analysis (RNA vs ATAC)
# ==============================================================================

cat("[Step 4] Computing cell type correlations...\n")

# --- RNA Correlation ---
DefaultAssay(seu) <- "RNA"
if(all(GetAssayData(seu, "RNA", "data") == 0)) seu <- NormalizeData(seu)

rna_data <- GetAssayData(seu, "RNA", "data")
# Select genes expressed in at least 1% of cells
expressed_genes <- rownames(rna_data)[rowSums(rna_data > 0) >= ncol(seu) * 0.01]

# Calculate average expression per cell type
rna_avg <- sapply(ordered_celltypes, function(ct) {
  cells <- colnames(seu)[seu$celltype == ct]
  if(length(cells) > 0) rowMeans(rna_data[expressed_genes, cells, drop=FALSE]) else NA
})
rna_cor <- cor(rna_avg, method = "spearman", use = "complete.obs")

# --- ATAC Correlation ---
DefaultAssay(seu) <- "ATAC"
atac_counts <- GetAssayData(seu, "ATAC", "counts")
# Select top 10k variable peaks
peak_var <- apply(atac_counts, 1, var)
top_peaks <- head(names(sort(peak_var, decreasing = TRUE)), 10000)

# Calculate average accessibility per cell type
atac_avg <- sapply(ordered_celltypes, function(ct) {
  cells <- colnames(seu)[seu$celltype == ct]
  if(length(cells) > 0) rowMeans(atac_counts[top_peaks, cells, drop=FALSE]) else NA
})
atac_cor <- cor(atac_avg, method = "spearman", use = "complete.obs")

# --- Combine Correlations ---
# Upper triangle: ATAC, Lower triangle: RNA, Diagonal: 1
combined_cor <- matrix(NA, nrow = n_celltypes, ncol = n_celltypes, dimnames = list(ordered_celltypes, ordered_celltypes))
for(i in 1:n_celltypes) {
  for(j in 1:n_celltypes) {
    if(i == j) combined_cor[i, j] <- 1
    else if(i < j) combined_cor[i, j] <- rna_cor[i, j] # Keep consistent with logic
    else combined_cor[i, j] <- atac_cor[i, j]
  }
}
# Note: Logic in original script put RNA in Lower (i<j depends on loop structure), 
# enforcing Upper=ATAC, Lower=RNA for clarity in heatmap.

# --- Plot Correlation Heatmap ---
pdf(file.path(OUTPUT_BASE, "Correlation_Analysis/correlation_heatmap.pdf"), 
    width = n_celltypes * 0.4 + 5, height = n_celltypes * 0.4 + 1)

ht_cor <- Heatmap(
  combined_cor,
  name = "Correlation",
  col = get_correlation_colors(),
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = FALSE, show_column_names = FALSE,
  border = TRUE,
  # Annotations
  top_annotation = HeatmapAnnotation(
    CellType = ordered_celltypes, col = list(CellType = color_dict),
    show_legend = TRUE, annotation_name_side = "left"
  ),
  left_annotation = rowAnnotation(
    CellType = ordered_celltypes, col = list(CellType = color_dict),
    show_legend = FALSE
  ),
  width = unit(n_celltypes * 0.4, "inch"),
  height = unit(n_celltypes * 0.4, "inch")
)
draw(ht_cor, heatmap_legend_side = "right")
dev.off()
cat("   -> Correlation heatmap saved.\n")

# ==============================================================================
# 6. Differential Expression (DEG) & Visualization
# ==============================================================================

cat("[Step 5] Analyzing Differential Expression (RNA)...\n")

markers_path <- file.path(OUTPUT_BASE, "Data/all_celltype_markers.csv")
if(file.exists(markers_path)) {
  all_markers <- read.csv(markers_path, stringsAsFactors = FALSE)
} else {
  DefaultAssay(seu) <- "RNA"
  Idents(seu) <- "celltype"
  all_markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
  write.csv(all_markers, markers_path, row.names = FALSE)
}

# Select Top 3 genes per cell type
top_genes <- all_markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(n = 3, order_by = avg_log2FC) %>%
  pull(gene) %>% unique()

# Downsample cells for visualization (max 150 per type)
set.seed(123)
selected_cells <- unlist(lapply(ordered_celltypes, function(ct) {
  cells <- colnames(seu)[seu$celltype == ct]
  if(length(cells) > 150) sample(cells, 150) else cells
}))
cell_type_vec <- seu$celltype[selected_cells]

# --- Plot DEG Heatmap ---
deg_mat <- as.matrix(GetAssayData(seu, "RNA", "data")[top_genes, selected_cells])
deg_mat <- t(scale(t(deg_mat))) # Z-score
deg_mat[is.na(deg_mat)] <- 0

pdf(file.path(OUTPUT_BASE, "DEG_Analysis/DEGs_high_density_heatmap.pdf"), width = 14, height = 12)

ht_deg <- Heatmap(
  deg_mat,
  name = "Expression",
  col = get_expression_colors(),
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_split = factor(cell_type_vec, levels = ordered_celltypes),
  top_annotation = HeatmapAnnotation(CellType = cell_type_vec, col = list(CellType = color_dict), show_legend = TRUE),
  use_raster = TRUE
)
draw(ht_deg, heatmap_legend_side = "right")
dev.off()
cat("   -> DEG heatmap saved.\n")

# ==============================================================================
# 7. Differential Motif Analysis & Visualization
# ==============================================================================

cat("[Step 6] Analyzing Differential Motifs (chromVAR)...\n")

motif_markers_path <- file.path(OUTPUT_BASE, "Data/differential_motif_markers.csv")
if(file.exists(motif_markers_path)) {
  motif_markers <- read.csv(motif_markers_path, stringsAsFactors = FALSE)
} else {
  DefaultAssay(seu) <- "chromvar"
  Idents(seu) <- "celltype"
  motif_markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
  write.csv(motif_markers, motif_markers_path, row.names = FALSE)
}

# Select Top 3 motifs per cell type
top_motifs <- motif_markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(n = 3, order_by = avg_log2FC) %>%
  pull(gene) %>% unique()

if(length(top_motifs) > 0) {
  # Clean motif names for display
  clean_names <- top_motifs
  try({
    clean_names <- ConvertMotifID(seu, id = top_motifs)
    clean_names <- gsub("MA[0-9]+\\.[0-9]+_", "", clean_names) # Remove JASPAR ID prefix if desired
    clean_names <- tools::toTitleCase(tolower(clean_names))
  }, silent = TRUE)
  
  # Prepare Matrix
  motif_mat <- as.matrix(GetAssayData(seu, "chromvar", "data")[top_motifs, selected_cells])
  motif_mat <- t(scale(t(motif_mat)))
  motif_mat[is.na(motif_mat)] <- 0
  rownames(motif_mat) <- clean_names
  
  # --- Plot Motif Heatmap ---
  pdf(file.path(OUTPUT_BASE, "TF_Motif_Analysis/motif_activity_heatmap.pdf"), width = 14, height = 12)
  
  ht_motif <- Heatmap(
    motif_mat,
    name = "TF Score",
    col = get_tf_colors(),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8, fontface = "bold"),
    column_split = factor(cell_type_vec, levels = ordered_celltypes),
    top_annotation = HeatmapAnnotation(CellType = cell_type_vec, col = list(CellType = color_dict), show_legend = TRUE),
    use_raster = TRUE
  )
  draw(ht_motif, heatmap_legend_side = "right")
  dev.off()
  cat("   -> Motif heatmap saved.\n")
} else {
  cat("   -> No significant motifs found.\n")
}

# ==============================================================================
# 8. Save Final Data
# ==============================================================================

cat("[Step 7] Saving final datasets...\n")
write.csv(combined_cor, file.path(OUTPUT_BASE, "Data/combined_celltype_correlation.csv"))
saveRDS(seu, file.path(OUTPUT_BASE, "Data/brain_multiome_final.rds"))

cat("\n=== Analysis Complete ===\n")
cat(sprintf("Results directory: %s\n", OUTPUT_BASE))

# Print session info for reproducibility
cat("\n=== Session Info ===\n")
sessionInfo()
