#!/usr/bin/env Rscript

# ======== 1. CONFIGURATION ========
WORKING_DIR  <- "."
DATA_DIR     <- "results/08_FigR_Analysis/Data"
OUTPUT_DIR   <- "results/08_FigR_Analysis/Plots"
SEURAT_FILE  <- "data/brain_multiome_cleaned.rds"

# List of genes to visualize
TARGET_GENES <- c("Bcl11b", "Gab1", "Sox9", "Olig1", "Gpr183") 

# ======== 2. SETUP & HELPERS ========
setwd(WORKING_DIR)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "Genome_Tracks"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "Triple_UMAPs"), showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(viridis)
  library(scales)
  library(GenomicRanges)
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

# --- Helper: Seurat v5 Safe Data Access ---
get_expr <- function(obj, assay, gene, layer="data") {
  DefaultAssay(obj) <- assay
  tryCatch({
    val <- GetAssayData(obj, assay=assay, layer=layer)[gene, ]
  }, error = function(e) {
    # Fallback for v4/v5 differences
    val <- LayerData(obj, assay=assay, layer=layer)[gene, ]
  })
  return(as.numeric(val))
}

# --- Helper: Plot UMAP with Custom Data ---
plot_custom_umap <- function(coords, expr, title, palette="viridis") {
  df <- data.frame(UMAP1=coords[,1], UMAP2=coords[,2], val=expr)
  df <- df[order(df$val), ] # Plot high values on top
  
  cols <- if(palette=="purple") c("#F7F7F7", "#3F007D") else 
          if(palette=="heat") c("#FFFFFF", "#CB181D") else 
          if(palette=="green") c("#F7FCF5", "#00441B") else viridis(100)
          
  ggplot(df, aes(x=UMAP1, y=UMAP2, color=val)) +
    geom_point(size=0.2, alpha=0.8) +
    scale_color_gradientn(colors=cols, oob=scales::squish) +
    theme_void() + ggtitle(title) + theme(legend.position="bottom")
}

# ======== 3. LOAD DATA ========
cat("[1/3] Loading Analysis Results...\n")
if (!file.exists(file.path(DATA_DIR, "figR_results.rds"))) stop("Run Step 1 first.")

seu <- readRDS(SEURAT_FILE)
cisCorr <- read.csv(file.path(DATA_DIR, "cisCorr_raw.csv"))
dorcMat <- readRDS(file.path(DATA_DIR, "dorcMat_smoothed.rds"))

# Filter for significant links
cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)

# Calculate Gene Activity if missing (needed for Triple Plot)
if(!"aRNA" %in% names(seu@assays)) {
  cat("   - Calculating Gene Activity (aRNA)...\n")
  act <- GeneActivity(seu)
  seu[["aRNA"]] <- CreateAssayObject(counts = act)
  seu <- NormalizeData(seu, assay = "aRNA")
}

# ======== 4. TRIPLE UMAP VISUALIZATION ========
cat("[2/3] Generating Triple Comparison Plots (RNA/Activity/DORC)...\n")

# Get UMAP Coords
umap <- Embeddings(seu, reduction = ifelse("atac.umap" %in% names(seu@reductions), "atac.umap", "umap"))

for(gene in TARGET_GENES) {
  # Check availability
  has_rna <- gene %in% rownames(seu[["RNA"]])
  has_dorc <- gene %in% rownames(dorcMat)
  
  if(has_rna && has_dorc) {
    p1 <- plot_custom_umap(umap, get_expr(seu, "RNA", gene), paste(gene, "RNA"), "purple")
    p2 <- plot_custom_umap(umap, get_expr(seu, "aRNA", gene), paste(gene, "Activity"), "green")
    p3 <- plot_custom_umap(umap, dorcMat[gene, rownames(umap)], paste(gene, "DORC"), "heat")
    
    combined <- p1 | p2 | p3
    ggsave(file.path(OUTPUT_DIR, "Triple_UMAPs", paste0(gene, "_triple.pdf")), combined, width=15, height=5)
    cat(sprintf("   - Saved Triple Plot for %s\n", gene))
  }
}

# ======== 5. GENOME TRACK PLOTS (FIXED ARCS) ========
cat("[3/3] Generating Genome Track Plots...\n")

# Note: Ideally, you would source the complex plotting function from a utils file
# For this script, we assume the `plot_DORC_regulation_fixed` function is defined here or sourced.
# (Due to length, I will provide the usage logic. You can paste the function from your previous Script 3 here)

# [Insert function `plot_DORC_regulation_with_all_genes_fixed` here]
# [Insert function `parse_peak_coordinates` here]

# --- Simplified Execution Loop ---
for(gene in TARGET_GENES) {
  # Ensure the gene has links
  if(gene %in% cisCorr.filt$Gene) {
    tryCatch({
      # Call the complex plotting function (assumed loaded)
      # You need to paste the `plot_DORC_regulation_with_all_genes_fixed` function definition above this loop
      # in the actual file.
      
      # For brevity in this response, I am showing the call structure:
      # p <- plot_DORC_regulation_with_all_genes_fixed(
      #   gene_name = gene,
      #   cisCorr_data = cisCorr.filt,
      #   seurat_obj = seu,
      #   extend_bp = 50000
      # )
      
      # if(!is.null(p)) {
      #   ggsave(file.path(OUTPUT_DIR, "Genome_Tracks", paste0(gene, "_track.pdf")), p, width=12, height=8)
      #   cat(sprintf("   - Saved Genome Track for %s\n", gene))
      # }
      
      cat(sprintf("   - (Placeholder) Genome Track logic for %s would run here.\n", gene))
      
    }, error = function(e) cat(sprintf("   ! Error plotting %s: %s\n", gene, e$message)))
  }
}

cat("\n=== Visualization Complete ===\n")
