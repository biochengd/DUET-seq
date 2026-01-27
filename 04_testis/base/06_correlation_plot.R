# ==============================================================================
# DUET-SEQ Analysis Pipeline: Gene-Peak Correlation Visualization
# ==============================================================================
# Description: 
# This script generates three types of plots for specific marker genes:
# 1. Boxplot of Gene Expression across cell types.
# 2. Boxplot of Chromatin Accessibility (best linked peak) across cell types.
# 3. Scatter plot showing the correlation between Gene Expression and Accessibility.
#
# Steps:
# 1. Load pre-calculated smoothed expression/accessibility matrices and links.
# 2. Define plotting function and parameters (colors, themes).
# 3. Generate plots for selected marker genes.
# ==============================================================================

# 1. Library Loading
# ==============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(qs)
  library(glue)
})

# 2. Directory Setup & Data Loading
# ==============================================================================
# Set working directory relative to project root
# setwd("./src/rna/seurat_signac_pipeline2/01_overall") 

data_dir <- "./data_overall"
outs_dir <- "./outs_overall/correlation_plots"

if (!dir.exists(outs_dir)) dir.create(outs_dir, recursive = TRUE)

message("Loading processed data matrices and links...")
# Load data processed in previous steps (e.g., from 04_visualization.R)
means <- qread(file = glue("{data_dir}/means.qs"))           # Smoothed RNA matrix
means_peak <- qread(file = glue("{data_dir}/means_peak.qs")) # Smoothed ATAC matrix
final_links_sorted <- qread(file = glue("{data_dir}/final_links_sorted.qs")) # Best links

# Reconstruct annotation dataframe (assuming uniform window size used in smoothing)
# Note: Ensure these levels match your data generation step
celltype_levels <- c(
  "Spermatogonia_SPG", "Spermatocytes_SPC", "Round_spermatids", 
  "Elongating_spermatids", "Sertoli_cells", "Leydig_cells", 
  "Stromal_cells", "Peritubular_myoid_cells", "Macrophage", "Endothelial_cells"
)

# Infer cell types from column names of the means matrix if they follow a pattern,
# or reconstruct based on known window size (n=10) and original cell counts.
# For simplicity here, we assume 'means' columns are named or ordered correctly.
# If you saved 'anno_col2' in the previous step, load it here. 
# Otherwise, we reconstruct a placeholder or simple extraction:
# Extract cell type from column names (assuming format "ID_Celltype")
celltype_major <- gsub("^[0-9]+", "", colnames(means)) 
# Or if you saved the annotation object:
# anno_col2 <- qread(glue("{data_dir}/anno_col2.qs")) 
# Using a reconstruction based on the script logic:
anno_col2 <- data.frame(celltype_major = factor(celltype_major, levels = celltype_levels))


# 3. Visualization Configuration
# ==============================================================================

# Cell type colors
cell_type_colors <- c(
  "Spermatogonia_SPG" = '#FB8D3C', "Spermatocytes_SPC" = '#d6604d',
  "Round_spermatids" = "#D082AF", "Elongating_spermatids" = "#7B599C", 
  "Sertoli_cells" = "#00ACB3", "Leydig_cells" = "#65AB53",
  "Stromal_cells" = "#4667A8", "Peritubular_myoid_cells" = "#9CC0DD", 
  "Macrophage" = "#85BFB8", "Endothelial_cells" = '#820610'
)

# Standard ggplot theme
plot_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks.x = element_blank(),
    legend.position = 'none',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# 4. Plotting Function
# ==============================================================================

#' Generate and save gene-peak correlation plots
#' @param gene_name Gene symbol to plot
#' @param links_data Dataframe containing gene-peak links (must contain 'gene', 'peak', 'score')
#' @param expr_data Smoothed expression matrix (rows=genes, cols=meta-cells)
#' @param atac_data Smoothed accessibility matrix (rows=peaks, cols=meta-cells)
#' @param cell_anno Dataframe with 'celltype_major' column corresponding to matrix columns
#' @param colors_vec Named vector of colors for cell types
#' @param output_dir Directory to save plots
generate_gene_plots <- function(gene_name, links_data, expr_data, atac_data, cell_anno, colors_vec, output_dir) {
  
  # Find best peak for the gene
  gene_link <- links_data %>%
    filter(gene == gene_name) %>%
    arrange(desc(score)) %>%
    head(1)
  
  if (nrow(gene_link) == 0) {
    message(glue("Skipping {gene_name}: No link found."))
    return(NULL)
  }
  
  top_peak <- gene_link$peak
  
  # 1. Gene Expression Boxplot
  df_expr <- data.frame(value = expr_data[gene_name, ], celltype = cell_anno$celltype_major)
  p1 <- ggplot(df_expr, aes(x = celltype, y = value, fill = celltype)) +
    geom_boxplot(outlier.size = 0.1) +
    scale_fill_manual(values = colors_vec) +
    labs(title = glue("{gene_name} expression"), x = "", y = "Relative expression") +
    plot_theme
  
  ggsave(glue("{output_dir}/{gene_name}_expression.pdf"), p1, width = 6, height = 4)
  
  # 2. Chromatin Accessibility Boxplot
  df_atac <- data.frame(value = atac_data[top_peak, ], celltype = cell_anno$celltype_major)
  p2 <- ggplot(df_atac, aes(x = celltype, y = value, fill = celltype)) +
    geom_boxplot(outlier.size = 0.1) +
    scale_fill_manual(values = colors_vec) +
    labs(title = glue("{top_peak} accessibility"), x = "", y = "Chromatin accessibility") +
    plot_theme
  
  ggsave(glue("{output_dir}/{gene_name}_accessibility.pdf"), p2, width = 6, height = 4)
  
  # 3. Correlation Scatter Plot
  df_corr <- data.frame(gene_expr = expr_data[gene_name, ], atac = atac_data[top_peak, ])
  p3 <- ggplot(df_corr, aes(x = gene_expr, y = atac)) +
    geom_point(size = 0.7, alpha = 0.6) +
    geom_smooth(color = "grey", method = "lm", se = FALSE) +
    labs(
      title = glue("{gene_name} & {top_peak}"),
      x = "Gene expression",
      y = "Chromatin accessibility"
    ) +
    plot_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
    
  ggsave(glue("{output_dir}/{gene_name}_correlation.pdf"), p3, width = 4, height = 4)
  
  message(glue("Generated plots for {gene_name}."))
}

# 5. Execution
# ==============================================================================
message("Generating plots for selected marker genes...")

# Define list of genes to plot (can be expanded)
# Specific markers of interest from the analysis
marker_genes_to_plot <- c(
    "Cd63", "Gopc", "Sun5", "Spag6", "Tnp1", "Prm1", "Prm2", "Gapdhs", "Adam3",
    "Uchl1", "Stra8", "Sycp2", "Piwil1", "Acrv1", "Spaca1", "Sox9", "Clu"
)

# Loop through genes
for (gene in marker_genes_to_plot) {
  # Check if gene exists in data before plotting
  if (gene %in% rownames(means) && any(final_links_sorted$gene == gene)) {
      generate_gene_plots(
        gene_name    = gene,
        links_data   = final_links_sorted,
        expr_data    = means,
        atac_data    = means_peak,
        cell_anno    = anno_col2,
        colors_vec   = cell_type_colors,
        output_dir   = outs_dir
      )
  } else {
      message(glue("Skipping {gene}: Not found in expression matrix or links."))
  }
}

message("Visualization pipeline completed.")
