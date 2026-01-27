# ==============================================================================
# DUET-SEQ Analysis Pipeline: Pseudobulk Differential Expression Analysis
# ==============================================================================
# Description: 
# This script performs pseudobulk differential expression analysis between 
# "Early" (P7-P21) and "Late" (P28-10W) stages. It includes:
# 1. Differential Expression Analysis (Volcano Plot).
# 2. Gene Set Enrichment Analysis (GSEA) using GO Biological Processes.
# 3. Visualization of enriched pathways.
#
# Steps:
# 1. Load the processed Seurat object.
# 2. Define developmental stages (Early vs Late).
# 3. Perform DE analysis and generate a volcano plot.
# 4. Run GSEA and visualize top enriched pathways.
# ==============================================================================

# 1. Library Loading
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(qs)
  library(glue)
  library(clusterProfiler)
  library(tibble)
  library(ggrepel)
})

# 2. Directory Setup & Data Loading
# ==============================================================================
# Set working directory relative to project root
# setwd("./src/rna/seurat_signac_pipeline2/04_re/01_preprocessing_qc_clustering")

data_dir <- "../data"
outs_dir <- "../outs"

if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(outs_dir)) dir.create(outs_dir, recursive = TRUE)

message("Loading Seurat object...")
# Assuming 'seu_re.qs' is the object to be analyzed
seu_re <- qread(file = glue("{data_dir}/seu_re.qs"))

# Define Stage based on Group (Timepoint)
early_stage <- c("P7", "P14", "P21")
late_stage_grps <- c("P28", "P35", "10W")

seu_re$stage_timepoint <- ifelse(seu_re$group %in% early_stage, "early", "late")
seu_re$stage_timepoint <- factor(seu_re$stage_timepoint, levels = c("early", "late"))

# Save updated metadata
qsave(seu_re, file = glue("{data_dir}/seu_re_processed.qs"))

# 3. Differential Expression Analysis (Volcano Plot)
# ==============================================================================
message("Performing Differential Expression Analysis (Early vs Late)...")

source("R/01_pseudobulk.R") # Assuming this helper script exists for plotting

DefaultAssay(seu_re) <- "SCT"
Idents(seu_re) <- "stage_timepoint"

# Find Markers
deg_results <- FindMarkers(
  seu_re,
  ident.1 = "late",       # Case
  ident.2 = "early",      # Control
  logfc.threshold = 0,    # Keep all for GSEA
  min.pct = 0,
  test.use = "wilcox",
  verbose = FALSE
) %>% rownames_to_column("gene")

# Filter noise genes
noise_pattern <- "^ENSMUS|^[0-9]|^Gm[0-9]|^Rik$|Rik[0-9]|Rik$|^LOC|^BC[0-9]|^mt-|^AC[0-9]+"
deg_clean <- deg_results[!grepl(noise_pattern, deg_results$gene), ]

# Plot Hyperbolic Volcano (Custom Function assumed from source)
# If function is not available, standard volcano plot logic can be used here.
tryCatch({
    p_volcano <- plot_hyperbolic_volcano_improved(deg_clean, logFC_threshold = 1.5, label_top_n = 25)
    ggsave(filename = glue("{outs_dir}/volcano_early_vs_late.pdf"), plot = p_volcano, width = 8, height = 8)
}, error = function(e) { message("Volcano plot function not found or failed.") })

# 4. Gene Set Enrichment Analysis (GSEA)
# ==============================================================================
message("Running GSEA...")

# Prepare ranked gene list
deg_sorted <- deg_clean %>% arrange(desc(avg_log2FC), p_val)
gene_list <- deg_sorted$avg_log2FC
names(gene_list) <- deg_sorted$gene
gene_list <- na.omit(gene_list)

# Load GMT file
gmt_file <- "./resource/m5.go.bp.v2025.1.Mm.symbols.gmt" # Update path as needed
if (file.exists(gmt_file)) {
    pathways <- read.gmt(gmt_file)
    
    # Run GSEA
    gsea_res <- GSEA(
        gene_list,
        TERM2GENE = pathways,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        minGSSize = 10,
        maxGSSize = 500,
        BPPARAM = BiocParallel::SerialParam()
    )
    
    # Save results
    result_df <- as.data.frame(gsea_res)
    qsave(result_df, file = glue("{data_dir}/gsea_results.qs"))
    
    # 5. Visualization of Enriched Pathways
    # ==============================================================================
    message("Visualizing GSEA results...")
    
    result_df <- result_df %>% mutate(log10_P = -log10(pvalue))
    
    # Top 10 Activated (Late)
    top_late <- result_df %>% 
        filter(NES > 0, p.adjust < 0.05) %>% 
        arrange(pvalue) %>% 
        head(10)
    
    # Top 10 Suppressed (Early)
    top_early <- result_df %>% 
        filter(NES < 0, p.adjust < 0.05) %>% 
        arrange(pvalue) %>% 
        head(10)
    
    # Plot Early Enriched (Suppressed in Late)
    color_early <- "#2E8B57" # Green
    p_early <- ggplot(top_early, aes(x = log10_P, y = reorder(Description, log10_P))) +
      geom_bar(stat = "identity", fill = color_early, width = 0.7) +
      scale_x_reverse(position = "top") +
      scale_y_discrete(position = "left") +
      labs(x = bquote(-log[10](P)), y = NULL, title = "Enriched in Early Stage") +
      theme_classic() +
      theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 10)
      )
    
    ggsave(filename = glue("{outs_dir}/gsea_early_enriched.pdf"), plot = p_early, width = 9, height = 9)
    
    # Plot Late Enriched (Activated in Late)
    color_late <- "#FF8C00"  # Orange
    p_late <- ggplot(top_late, aes(x = log10_P, y = reorder(Description, log10_P))) +
      geom_bar(stat = "identity", fill = color_late, width = 0.7) +
      scale_x_continuous(position = "top") +
      scale_y_discrete(position = "right") +
      labs(x = bquote(-log[10](P)), y = NULL, title = "Enriched in Late Stage") +
      theme_classic() +
      theme(
        axis.text.y = element_text(hjust = 0, size = 10),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    
    ggsave(filename = glue("{outs_dir}/gsea_late_enriched.pdf"), plot = p_late, width = 9, height = 9)
    
    # Combined Plot
    p_combined <- p_late | p_early
    ggsave(filename = glue("{outs_dir}/gsea_combined.pdf"), plot = p_combined, width = 18, height = 9)
    
} else {
    message("GMT file not found. Skipping GSEA.")
}

message("Pseudobulk analysis completed.")
