# --- 1. Load Libraries ---
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(scales)

set.seed(1234)

# --- 2. Global Settings & Helper Functions ---
output_dir <- "/s1/chengd/project/multiome/02_cellline/cellline_trackplot"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Helper to filter standard chromosomes
Standard_peaks <- function(rowRanges, genome_type = 'hg38') {
  if (genome_type == 'hg38') {
    main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
  } else {
    main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
  }
  keep.peaks <- which(as.character(seqnames(rowRanges)) %in% main.chroms)
  return(rowRanges[keep.peaks])
}

# Helper to process BED files into GRanges
Processing_peaks <- function(peak_file, genome_type = 'hg38') {
  peaks_df <- read.table(peak_file, col.names = c("chr", "start", "end"))
  gr.peaks <- makeGRangesFromDataFrame(peaks_df)
  gr.peaks <- Standard_peaks(gr.peaks, genome_type = genome_type)
  return(gr.peaks)
}

# --- 3. Extract Gene Annotations ---
message(">>> Extracting Human (hg38) and Mouse (mm10) annotations...")
# Human
annotations_hg38 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations_hg38) <- 'UCSC'
genome(annotations_hg38) <- "hg38"

# Mouse
annotations_mm10 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations_mm10) <- 'UCSC'
genome(annotations_mm10) <- "mm10"


# ==============================================================================
# SECTION A: Human 293T Comparison
# ==============================================================================
message(">>> Starting Human 293T Processing...")

# --- A.1 Path Configuration ---
duet_293t_base <- "/s1/chengd/project/multiome/01_preprocessing/duet-seq/293t/293T_chromap/outs/"
issaac_293t_base <- "/s1/chengd/project/multiome/01_preprocessing/issaac-seq/ISSAAC_ATAC_Analysis/outs/"

# --- A.2 Object Construction ---
peaks_293t <- Processing_peaks(paste0(duet_293t_base, "raw_mtx/peaks.bed"), 'hg38')

# Duet-seq 293T
duet_cells_293t <- readLines(paste0(duet_293t_base, "filtered_mtx/barcodes.tsv"))
duet_frag_293t <- CreateFragmentObject(path = paste0(duet_293t_base, "fragments.tsv.gz"), cells = duet_cells_293t)
duet_counts_293t <- FeatureMatrix(fragments = duet_frag_293t, features = peaks_293t, cells = duet_cells_293t)
duet_293t_obj <- CreateSeuratObject(counts = duet_counts_293t, assay = "ATAC")
duet_293t_obj[["ATAC"]] <- CreateChromatinAssay(counts = duet_counts_293t, fragments = duet_frag_293t, genome = "hg38")
Annotation(duet_293t_obj) <- annotations_hg38

# ISSAAC-seq 293T (Handling GRCh38 prefix)
issaac_cells_293t <- readLines(paste0(issaac_293t_base, "filtered_mtx/barcodes.tsv"))
issaac_cells_293t <- gsub("-1$", "", issaac_cells_293t[issaac_cells_293t != ""])
issaac_frag_293t <- CreateFragmentObject(path = paste0(issaac_293t_base, "fragments.tsv.gz"), cells = issaac_cells_293t)

peaks_issaac_prefix <- renameSeqlevels(peaks_293t, paste0("GRCh38_", seqlevels(peaks_293t)))
issaac_counts_293t <- FeatureMatrix(fragments = issaac_frag_293t, features = peaks_issaac_prefix, cells = issaac_cells_293t)
rownames(issaac_counts_293t) <- gsub("^GRCh38_", "", rownames(issaac_counts_293t))

issaac_293t_obj <- CreateSeuratObject(counts = issaac_counts_293t, assay = "ATAC")
issaac_293t_obj[["ATAC"]] <- CreateChromatinAssay(counts = issaac_counts_293t, fragments = issaac_frag_293t, genome = "hg38")
Annotation(issaac_293t_obj) <- annotations_hg38

# --- A.3 Visualization ---
region_293t <- "chr11-65600000-66000000"
region_293t_prefix <- paste0("GRCh38_", region_293t)

p1 <- CoveragePlot(duet_293t_obj, region = region_293t, annotation = F, peaks = F) + scale_fill_manual(values = "#E41A1C") + labs(y = "Duet ATAC")
p2 <- CoveragePlot(issaac_293t_obj, region = region_293t_prefix, annotation = F, peaks = F) + scale_fill_manual(values = "#4DAF4A") + labs(y = "ISSAAC ATAC")
p3 <- TilePlot(duet_293t_obj, region = region_293t, tile.cells = 300) + 
      scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1), oob = squish) + theme(legend.position = "none") + labs(y = "Duet Tile")
p4 <- TilePlot(issaac_293t_obj, region = region_293t_prefix, tile.cells = 200) + 
      scale_fill_gradient(low = "white", high = "#4DAF4A", limits = c(0, 1), oob = squish) + theme(legend.position = "none") + labs(y = "ISSAAC Tile")
p5 <- AnnotationPlot(duet_293t_obj, region = region_293t)

combined_293t <- (p1 / p2 / p3 / p4 / p5) + plot_layout(heights = c(2, 2, 3, 3, 1))
ggsave(file.path(output_dir, "293T_Duet_vs_ISSAAC.pdf"), combined_293t, width = 12, height = 16)


# ==============================================================================
# SECTION B: Mouse 3T3 Comparison
# ==============================================================================
message(">>> Starting Mouse 3T3 Processing...")

# --- B.1 Path Configuration ---
duet_3t3_base <- "/s1/chengd/project/multiome/01_preprocessing/duet-seq/3t3/3T3_chromap/outs/"
tenx_3t3_base <- "/s1/chengd/project/multiome/01_preprocessing/10x_multiome/3t3/NIH3T3_ATAC_Analysis/outs/"
rna_bw_duet <- "/s1/chengd/project/multiome/01_preprocessing/duet-seq/3t3/3T3_rna/outs/rna_coverage.bw"
rna_bw_10x <- "/s1/chengd/project/multiome/01_preprocessing/10x_multiome/3t3/rna/10x_rna_coverage.bw"

# --- B.2 Object Construction ---
peaks_3t3 <- Processing_peaks(paste0(duet_3t3_base, "raw_mtx/peaks.bed"), 'mm10')
peaks_3t3 <- peaks_3t3[width(peaks_3t3) < 5000 & width(peaks_3t3) > 100]

# Duet-seq 3T3
duet_cells_3t3 <- readLines(paste0(duet_3t3_base, "filtered_mtx/barcodes.tsv"))
duet_frag_3t3 <- CreateFragmentObject(path = paste0(duet_3t3_base, "fragments.tsv.gz"), cells = duet_cells_3t3)
duet_counts_3t3 <- FeatureMatrix(fragments = duet_frag_3t3, features = peaks_3t3, cells = duet_cells_3t3)
duet_3t3_obj <- CreateSeuratObject(counts = duet_counts_3t3, assay = "ATAC")
duet_3t3_obj[["ATAC"]] <- CreateChromatinAssay(counts = duet_counts_3t3, fragments = duet_frag_3t3, genome = "mm10")
Annotation(duet_3t3_obj) <- annotations_mm10

# 10x Multiome 3T3
tenx_cells_3t3 <- readLines(paste0(tenx_3t3_base, "filtered_mtx/barcodes.tsv"))
tenx_frag_3t3 <- CreateFragmentObject(path = paste0(tenx_3t3_base, "fragments.tsv.gz"), cells = tenx_cells_3t3)
tenx_counts_3t3 <- FeatureMatrix(fragments = tenx_frag_3t3, features = peaks_3t3, cells = tenx_cells_3t3)
tenx_3t3_obj <- CreateSeuratObject(counts = tenx_counts_3t3, assay = "ATAC")
tenx_3t3_obj[["ATAC"]] <- CreateChromatinAssay(counts = tenx_counts_3t3, fragments = tenx_frag_3t3, genome = "mm10")
Annotation(tenx_3t3_obj) <- annotations_mm10

# --- B.3 Visualization ---
region_3t3 <- "chr11-68900000-69000000"

m1 <- BigwigTrack(region = region_3t3, bigwig = rna_bw_duet, y_label = "Duet RNA") + scale_fill_manual(values = "#4682B4")
m2 <- BigwigTrack(region = region_3t3, bigwig = rna_bw_10x, y_label = "10x RNA") + scale_fill_manual(values = "#5F9EA0")
m3 <- CoveragePlot(duet_3t3_obj, region = region_3t3, annotation = F, peaks = F) + scale_fill_manual(values = "#E41A1C") + labs(y = "Duet ATAC")
m4 <- CoveragePlot(tenx_3t3_obj, region = region_3t3, annotation = F, peaks = F) + scale_fill_manual(values = "#377EB8") + labs(y = "10x ATAC")
m5 <- TilePlot(duet_3t3_obj, region = region_3t3, tile.cells = 200, order.by = "total") + 
      scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1), oob = squish) + 
      theme(legend.position = "none") + labs(y = "200 cells\n(Duet ATAC)")
m6 <- AnnotationPlot(duet_3t3_obj, region = region_3t3)

combined_3t3 <- (m1 / m2 / m3 / m4 / m5 / m6) + plot_layout(heights = c(1.5, 1.5, 2, 2, 4, 1))
ggsave(file.path(output_dir, "3T3_ATAC_RNA_Comparison.pdf"), combined_3t3, width = 12, height = 18)

message(">>> All tasks completed successfully.")
