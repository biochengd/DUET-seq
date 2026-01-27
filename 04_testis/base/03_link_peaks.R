# ==============================================================================
# DUET-SEQ Analysis Pipeline: Peak-Gene Linkage Analysis
# ==============================================================================
# Description: 
# This script identifies cis-regulatory elements by linking ATAC-seq peaks 
# to gene expression using the Signac LinkPeaks function. It focuses on 
# cell-type-specific marker genes identified in previous steps.
#
# Steps:
# 1. Load the fully processed and annotated Seurat object.
# 2. Identify cell-type specific marker genes.
# 3. Filter marker genes (remove noise/low specificity).
# 4. Calculate GC content and other region stats for peaks.
# 5. Link peaks to gene expression.
# ==============================================================================

# 1. Library Loading
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(qs)
  library(glue)
  library(dplyr)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(EnsDb.Mmusculus.v79)
  library(GenomicRanges)
})

# 2. Directory Setup & Data Loading
# ==============================================================================
# Set working directory relative to project root
# setwd("./src/rna/seurat_signac_pipeline2") 

data_dir <- "./data_overall"
outs_dir <- "./outs_overall"

if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(outs_dir)) dir.create(outs_dir, recursive = TRUE)

message("Loading annotated Seurat object...")
# Load the object from the previous step (e.g., after leiden clustering and annotation)
seu_obj <- qread(file = glue("{data_dir}/seu_processed_annotated.qs"))

# 3. Identify Cell-Type Specific Marker Genes
# ==============================================================================
message("Identifying cell-type specific marker genes...")

DefaultAssay(seu_obj) <- "SCT"
Idents(seu_obj) <- "celltype"

# Prepare SCT assay for differential expression (optional but recommended)
seu_obj <- PrepSCTFindMarkers(seu_obj)

# Find markers for all clusters
# Criteria: positive markers only, min 10% cells expressing, logFC threshold 0.25
de_genes <- FindAllMarkers(
  seu_obj, 
  only.pos = TRUE, 
  min.pct = 0.1, 
  logfc.threshold = 0.25,
  verbose = FALSE
)

# 4. Filter Marker Genes
# ==============================================================================
message("Filtering marker genes...")

# Initial filter: Adjusted P-value < 0.05 and Log2FC > 0.1
sig_markers <- de_genes %>% 
  filter(p_val_adj < 0.05, avg_log2FC > 0.1)

# Remove noise genes (Ensembl IDs, Rik genes, Gm genes, etc.)
noise_pattern <- "^ENSMUS|^[0-9]|^Gm[0-9]|^Rik$|Rik[0-9]|Rik$|^LOC|^BC[0-9]"
clean_markers <- sig_markers[!grepl(noise_pattern, sig_markers$gene), ]

# Specificity filter: Difference in percentage expression (pct.1 - pct.2) > 0.1
clean_markers$specificity <- clean_markers$pct.1 - clean_markers$pct.2
final_markers <- clean_markers %>% filter(specificity > 0.1)

# Extract unique gene list for linking
genes_to_link <- unique(final_markers$gene)

message(glue("Identified {length(genes_to_link)} specific marker genes for linkage analysis."))

# Save marker results
qsave(final_markers, file = glue("{data_dir}/final_celltype_markers.qs"))

# 5. Link Peaks to Genes
# ==============================================================================
message("Calculating peak-gene links (this may take time)...")

DefaultAssay(seu_obj) <- "ATAC"

# 5.1 Calculate Region Statistics (GC content)
# Required for background peak selection in LinkPeaks
seu_obj <- RegionStats(
  seu_obj, 
  assay = "ATAC", 
  genome = BSgenome.Mmusculus.UCSC.mm10
)

# 5.2 Run LinkPeaks
# Links peaks within 500kb of the TSS to gene expression
seu_obj <- LinkPeaks(
  object = seu_obj,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = genes_to_link,
  distance = 500000,   # 500 kb window
  min.cells = 10,      # Minimum cells expressing gene/peak
  n_sample = 200,      # Downsample for p-value calculation (speeds up process)
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
  verbose = TRUE
)

# 6. Save Results
# ==============================================================================
message("Saving results...")

# Save the updated Seurat object with links
qsave(seu_obj, file = glue("{data_dir}/seu_with_links.qs"))

# Extract and save the Links GRanges object separately for easier access
links_granges <- Links(seu_obj[["ATAC"]])
qsave(links_granges, file = glue("{data_dir}/peak_gene_links.qs"))

message("Linkage analysis completed successfully.")
