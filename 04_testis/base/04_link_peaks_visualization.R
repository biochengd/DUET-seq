# ==============================================================================
# DUET-SEQ Analysis Pipeline: Visualization (Heatmaps & Coverage Plots)
# ==============================================================================
# Description: 
# This script generates visualizations for the Multiome analysis, including:
# 1. Paired heatmaps for Gene Expression (RNA) and Chromatin Accessibility (ATAC).
# 2. Coverage plots for specific marker genes.
#
# Steps:
# 1. Load the processed Seurat object and LinkPeaks results.
# 2. Prepare data for heatmaps (smoothing, ordering).
# 3. Generate and save heatmaps.
# 4. Generate and save CoveragePlots for key genes.
# ==============================================================================

# 1. Library Loading
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(pheatmap)
  library(viridis)
  library(qs)
  library(glue)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(EnsDb.Mmusculus.v79)
  library(zoo) # For rolling mean smoothing
})

# 2. Directory Setup & Data Loading
# ==============================================================================
# Set working directory relative to project root
# setwd("./src/rna/seurat_signac_pipeline2") 

data_dir <- "./data_overall"
outs_dir <- "./outs_overall"

if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(outs_dir)) dir.create(outs_dir, recursive = TRUE)

message("Loading data...")
# Load Seurat object and LinkPeaks results from previous steps
seu_obj <- qread(file = glue("{data_dir}/seu_with_links.qs"))
links_granges <- qread(file = glue("{data_dir}/peak_gene_links.qs"))
# Load marker genes if needed, or re-calculate/extract from object
# Assuming markers were saved or can be extracted
# For this script, we'll re-use the logic to define specific genes/peaks

# Define standard cell type levels for ordering
celltype_levels <- c(
  "Spermatogonia_SPG", "Spermatocytes_SPC", "Round_spermatids", 
  "Elongating_spermatids", "Sertoli_cells", "Leydig_cells", 
  "Stromal_cells", "Peritubular_myoid_cells", "Macrophage", "Endothelial_cells"
)

# Ensure cell types in object match these levels
# Handle potential name changes (e.g., Telocyte -> Stromal_cells)
if ("Telocyte" %in% levels(seu_obj$celltype)) {
    seu_obj$celltype <- gsub("Telocyte", "Stromal_cells", seu_obj$celltype)
}
seu_obj$celltype <- factor(seu_obj$celltype, levels = celltype_levels)

# 3. Process Links for Heatmap
# ==============================================================================
message("Processing links for heatmap...")

links_df <- as.data.frame(mcols(links_granges))
links_df$adj_pval <- p.adjust(links_df$pvalue, method = 'BH')

# Filter links: Significant and Positive correlation
sig_links <- links_df %>% filter(adj_pval < 0.05, score > 0)

# Assign primary cell type to each gene based on marker specificity
# (This logic assumes you have a table 'top_markers_per_celltype' or similar)
# For the purpose of this clean script, we assume we want to plot links 
# for genes that are markers for the defined cell types.

# Identify markers (if not loaded) - simplified for script
DefaultAssay(seu_obj) <- "SCT"
Idents(seu_obj) <- "celltype"
markers <- FindAllMarkers(seu_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)

# Map genes to their top cell type
gene_to_celltype <- markers %>%
  group_by(gene) %>%
  slice_max(n = 1, order_by = avg_log2FC, with_ties = FALSE) %>%
  select(gene, primary_celltype = cluster)

# Join with links
links_annotated <- sig_links %>%
  inner_join(gene_to_celltype, by = "gene") %>%
  filter(!is.na(primary_celltype))

# Select best peak per gene for visualization
best_links <- links_annotated %>%
  group_by(gene) %>%
  slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(factor(primary_celltype, levels = celltype_levels), desc(zscore))

genes_to_plot <- best_links$gene
peaks_to_plot <- best_links$peak

message(glue("Selected {length(genes_to_plot)} gene-peak pairs for heatmap."))

# 4. Prepare Data for Heatmap (Smoothing)
# ==============================================================================
message("Smoothing data for heatmap...")

# Order cells by cell type
cell_order <- seu_obj@meta.data %>%
  arrange(celltype) %>%
  rownames()

# Helper function for smoothing
smooth_matrix <- function(mat, cell_groups, n = 10) {
  smoothed_list <- lapply(cell_groups, function(cells) {
    if (length(cells) >= n) {
      sub_mat <- mat[, cells, drop = FALSE]
      num_windows <- floor(length(cells) / n)
      # Calculate mean for each window of n cells
      binned_means <- t(apply(sub_mat[, 1:(num_windows * n)], 1, function(row) {
          sapply(1:num_windows, function(i) mean(row[((i-1)*n + 1):(i*n)]))
      }))
      return(binned_means)
    }
    return(NULL)
  })
  do.call(cbind, Filter(Negate(is.null), smoothed_list))
}

# Extract and smooth RNA
rna_data <- GetAssayData(seu_obj, assay = "SCT", layer = "scale.data")[genes_to_plot, cell_order]
cells_by_type <- split(cell_order, seu_obj$celltype[cell_order])
rna_smoothed <- smooth_matrix(rna_data, cells_by_type)

# Extract and smooth ATAC
atac_data <- GetAssayData(seu_obj, assay = "ATAC", layer = "counts")[peaks_to_plot, cell_order]
atac_norm <- RunTFIDF(atac_data) # Apply TF-IDF
atac_smoothed <- smooth_matrix(atac_norm, cells_by_type)

# Prepare Annotation
meta_cell_info <- do.call(rbind, lapply(names(cells_by_type), function(ct) {
    n_windows <- floor(length(cells_by_type[[ct]]) / 10)
    if (n_windows > 0) return(data.frame(celltype = rep(ct, n_windows)))
    return(NULL)
}))
rownames(meta_cell_info) <- colnames(rna_smoothed)

# 5. Generate Heatmaps
# ==============================================================================
message("Generating heatmaps...")

# Colors
cell_type_colors <- c(
  "Spermatogonia_SPG" = '#FB8D3C', "Spermatocytes_SPC" = '#d6604d',
  "Round_spermatids" = "#D082AF", "Elongating_spermatids" = "#7B599C", 
  "Sertoli_cells" = "#00ACB3", "Leydig_cells" = "#65AB53",
  "Stromal_cells" = "#4667A8", "Peritubular_myoid_cells" = "#9CC0DD", 
  "Macrophage" = "#85BFB8", "Endothelial_cells" = '#820610'
)
ann_colors <- list(celltype = cell_type_colors)

# Highlight specific genes
highlight_genes <- c("Sox9", "Gata4", "Stra8", "Sycp3", "Acrv1", "Prm1", "Dazl", "Utf1")
labels_row <- rownames(rna_smoothed)
labels_row[!labels_row %in% highlight_genes] <- ""

# RNA Heatmap
pdf(glue("{outs_dir}/Gene_Expression_Heatmap.pdf"), width = 12, height = 10, useDingbats = FALSE)
pheatmap(
  rna_smoothed,
  scale = "row",
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = TRUE, show_colnames = FALSE,
  labels_row = labels_row,
  annotation_col = meta_cell_info,
  annotation_colors = ann_colors,
  color = viridis(100),
  main = "Gene Expression",
  use_raster = TRUE
)
dev.off()

# ATAC Heatmap
pdf(glue("{outs_dir}/Chromatin_Accessibility_Heatmap.pdf"), width = 12, height = 10, useDingbats = FALSE)
pheatmap(
  atac_smoothed,
  scale = "row",
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  annotation_col = meta_cell_info,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#132B43", "#56B1F7", "#FFFFBF", "#FFED36", "#E31A1C"))(100),
  main = "Chromatin Accessibility",
  use_raster = TRUE
)
dev.off()
