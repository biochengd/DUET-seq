#!/usr/bin/env Rscript

# ======== 1. CONFIGURATION ========
# Working Directory
WORKING_DIR <- "/s1/chengd/project/multiome/03_brain/brain_multiome_v2_N1N2"

# Input Seurat Object (.rds)
INPUT_RDS <- "03_Cleaned_Analysis/Results/brain_multiome_cleaned.rds"

# Output Directory
OUTPUT_DIR <- "03_Cleaned_Analysis/Gene_Trackplots"

# List of genes to visualize
TARGET_GENES <- c("Cst3") 

# Plotting Parameters
EXTEND_BP <- 20000    # Flanking region size (bp)
TILE_CELLS <- 500     # Max cells to display in TilePlot

# ==============================================================================

# Setup Environment
setwd(WORKING_DIR)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(ggplot2)
  library(patchwork)
})

# ======== 2. DATA LOADING & SETUP ========
cat("Loading Seurat object...\n")
seu <- readRDS(INPUT_RDS)
DefaultAssay(seu) <- "ATAC"

# Setup Genome Annotations (Load once for efficiency)
cat("Setting up genome annotations...\n")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"
Annotation(seu) <- annotations

# ======== 3. PLOTTING FUNCTION ========
generate_gene_report <- function(obj, gene, out_dir, extend_bp, tile_cnt) {
  
  cat(sprintf("Processing gene: %s ... ", gene))
  
  tryCatch({
    # 1. Get gene coordinates
    region <- LookupGeneCoords(obj, gene = gene)
    if(is.null(region)) stop("Gene coordinates not found.")
    
    # Extend region
    region <- GRanges(
      seqnames = seqnames(region),
      ranges = IRanges(start = max(1, start(region) - extend_bp), 
                       end = end(region) + extend_bp)
    )
    
    # 2. Coverage Plot (Grouped by Cell Type)
    p_cov <- CoveragePlot(
      object = obj,
      region = region,
      group.by = "celltype",
      annotation = TRUE,
      peaks = TRUE,
      links = FALSE # Set to TRUE if LinkPeaks has been run
    ) + ggtitle(paste0(gene, " - Chromatin Accessibility by Cell Type"))
    
    # 3. Tile Plot
    p_tile <- TilePlot(
      object = obj,
      region = region,
      group.by = "celltype",
      tile.cells = tile_cnt
    ) + scale_fill_gradient(low = "white", high = "black") + 
      ggtitle("Fragment Density") + NoLegend()

    # 4. Combine and Save
    combined_plot <- p_cov / p_tile + plot_layout(heights = c(3, 1))
    
    out_file <- file.path(out_dir, paste0(gene, "_genomic_tracks.pdf"))
    ggsave(out_file, combined_plot, width = 12, height = 10)
    
    cat("Done. Saved to", out_file, "\n")
    
  }, error = function(e) {
    cat("Failed.\nError: ", e$message, "\n")
  })
}

# ======== 4. EXECUTION ========
cat("Starting batch processing for", length(TARGET_GENES), "genes...\n")

for (gene in TARGET_GENES) {
  generate_gene_report(seu, gene, OUTPUT_DIR, EXTEND_BP, TILE_CELLS)
}

cat("=== Analysis Complete ===\n")
