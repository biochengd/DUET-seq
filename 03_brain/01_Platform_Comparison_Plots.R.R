#!/usr/bin/env Rscript

# ======== 1. CONFIGURATION ========
WORKING_DIR  <- "/s1/chengd/project/multiome/03_brain/brain_multiome_v1_N1N210X"
INPUT_FILE   <- "03_Cleaned_Analysis/Results/brain_multiome_cleaned_final.rds"
OUTPUT_DIR   <- "04_Final_Figures"
N_CORES      <- 8

# ==============================================================================

setwd(WORKING_DIR)

# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(future)
})

# Parallel setup
plan("multisession", workers = N_CORES)
options(future.globals.maxSize = 80000 * 1024^2)

# Create specific directories
dir.create(OUTPUT_DIR, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "QC_comparison"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "UMAP_integration"), showWarnings = FALSE)

# ======== 2. LOAD DATA ========
cat("[1/5] Loading data...\n")
if (!file.exists(INPUT_FILE)) stop("Input file not found.")
seu <- readRDS(INPUT_FILE)
cat(sprintf("   - Loaded %d cells.\n", ncol(seu)))

# Helper: Add grouping column for replicates
seu$platform_replicate <- paste0(seu$platform, "_", seu$replicate)

# ======== 3. RNA QC METRICS ========
cat("[2/5] Generating RNA QC Plots...\n")

rna_metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
plots_rna <- list()

for (m in rna_metrics) {
  p <- VlnPlot(seu, features = m, group.by = "platform_replicate", pt.size = 0) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(x = NULL)
  
  if (m != "percent.mt") p <- p + scale_y_log10()
  plots_rna[[m]] <- p
}

# Clean Title: Removed "Figure A"
fig_rna <- wrap_plots(plots_rna, ncol = 3) + 
  plot_annotation(title = "RNA Quality Control Metrics", theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

# Clean Filename
ggsave(file.path(OUTPUT_DIR, "QC_comparison/RNA_QC_ViolinPlots.pdf"), fig_rna, width = 15, height = 5)

# ======== 4. ATAC QC METRICS ========
cat("[3/5] Generating ATAC QC Plots...\n")

atac_metrics <- c("nFeature_ATAC", "nCount_ATAC", "TSS.enrichment", "FRiP")
plots_atac <- list()

for (m in atac_metrics) {
  if(m %in% colnames(seu@meta.data)) {
    p <- VlnPlot(seu, features = m, group.by = "platform_replicate", pt.size = 0) +
      geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      labs(x = NULL)
    
    if (m %in% c("nFeature_ATAC", "nCount_ATAC")) p <- p + scale_y_log10()
    plots_atac[[m]] <- p
  }
}

# Clean Title: Removed "Figure B"
fig_atac <- wrap_plots(plots_atac, ncol = 4) + 
  plot_annotation(title = "ATAC Quality Control Metrics", theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

# Clean Filename
ggsave(file.path(OUTPUT_DIR, "QC_comparison/ATAC_QC_ViolinPlots.pdf"), fig_atac, width = 20, height = 5)

# ======== 5. COMBINED QC SUMMARY ========
cat("[4/5] Generating Combined QC Summary...\n")

# Prepare data manually for facet plotting
meta <- seu@meta.data
df_long <- rbind(
  data.frame(val = log10(meta$nFeature_RNA), metric = "Log10(Genes)", group = meta$platform_replicate, platform = meta$platform),
  data.frame(val = log10(meta$nCount_RNA), metric = "Log10(RNA UMIs)", group = meta$platform_replicate, platform = meta$platform),
  data.frame(val = log10(meta$nCount_ATAC), metric = "Log10(ATAC Reads)", group = meta$platform_replicate, platform = meta$platform)
)

p_comb <- ggplot(df_long, aes(x = group, y = val, fill = platform)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  facet_wrap(~metric, scales = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
  labs(x = NULL, y = "Value", fill = "Platform")

ggsave(file.path(OUTPUT_DIR, "QC_comparison/Combined_QC_Summary.pdf"), p_comb, width = 15, height = 6)

# ======== 6. UMAP INTEGRATION VISUALIZATION ========
cat("[5/5] Generating UMAP Plots...\n")

# Determine best reduction to plot (WNN > ATAC > RNA)
red_use <- if("wnn.umap" %in% names(seu@reductions)) "wnn.umap" else "umap"

# 1. Integration View (By Platform)
p_int <- DimPlot(seu, reduction = red_use, group.by = "platform", shuffle = TRUE, pt.size = 0.1) +
  ggtitle("Integration by Platform") + theme_void() + 
  theme(plot.title = element_text(hjust = 0.5))

# 2. Annotation View (By Cell Type)
p_anno <- DimPlot(seu, reduction = red_use, group.by = "celltype", label = TRUE, repel = TRUE, pt.size = 0.1) +
  ggtitle("Cell Type Annotation") + theme_void() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

# 3. Split View (Side-by-Side)
p_split <- DimPlot(seu, reduction = red_use, group.by = "celltype", split.by = "platform", pt.size = 0.1) +
  ggtitle("Split by Platform") + theme_void() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right")

# Save plots
ggsave(file.path(OUTPUT_DIR, "UMAP_integration/UMAP_Integration_ByPlatform.pdf"), p_int, width = 8, height = 7)
ggsave(file.path(OUTPUT_DIR, "UMAP_integration/UMAP_Annotation_ByCelltype.pdf"), p_anno, width = 8, height = 7)
ggsave(file.path(OUTPUT_DIR, "UMAP_integration/UMAP_Split_Comparison.pdf"), p_split, width = 14, height = 7)

cat("\n=== Analysis Complete ===\n")
cat("Generated Files in '04_Final_Figures/':\n")
cat("1. QC_comparison/RNA_QC_ViolinPlots.pdf\n")
cat("2. QC_comparison/ATAC_QC_ViolinPlots.pdf\n")
cat("3. QC_comparison/Combined_QC_Summary.pdf\n")
cat("4. UMAP_integration/ (3 files)\n")
