# ==============================================================================
# DUET-SEQ Analysis Pipeline: Spermatogonia (SPG) Sub-analysis
# ==============================================================================
# Description: 
# This script focuses on the Spermatogonia (SPG) population extracted from the 
# overall dataset. It performs:
# 1. Stage definition (Early vs Late).
# 2. Identification of marker genes for different timepoints and stages.
# 3. Visualization of SPG subtypes and developmental trajectories (Pseudotime).
# 4. Cellular composition analysis across pseudotime.
#
# Steps:
# 1. Load the SPG subset Seurat object.
# 2. Define developmental stages.
# 3. Find top marker genes for groups and stages.
# 4. Generate DotPlots for known markers.
# 5. Visualize Pseudotime distributions (RidgePlots).
# 6. Analyze cellular composition dynamics (Stacked Area Plot).
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
  library(ggrastr)
  library(RColorBrewer)
})

# 2. Directory Setup & Data Loading
# ==============================================================================
# Set working directory relative to project root
# setwd("./src/rna/seurat_signac_pipeline2/02_spg/01_preprocessing_qc_clustering")

data_dir <- "../data"
outs_dir <- "../outs"

if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(outs_dir)) dir.create(outs_dir, recursive = TRUE)

message("Loading SPG Seurat object...")
# Assuming 'seu_spg.qs' is the subsetted object containing only Spermatogonia
seu_spg <- qread(file = glue("{data_dir}/seu_spg.qs"))

# Define Stage based on Group (Timepoint)
# P7, P14 -> Early; Others -> Late
seu_spg$stage <- ifelse(seu_spg$group %in% c("P7", "P14"), "Early", "Late")

# Standardize group names (e.g., 10W -> W10) and set levels
seu_spg$group <- as.character(seu_spg$group)
seu_spg$group[seu_spg$group == "10W"] <- "W10"
seu_spg$group <- factor(seu_spg$group, levels = c("P7","P14", "P21", "P28", "P35", "W10"))

# Save updated metadata
qsave(seu_spg, file = glue("{data_dir}/seu_spg_processed.qs"))

# 3. Marker Gene Identification
# ==============================================================================
message("Identifying marker genes...")

DefaultAssay(seu_spg) <- "SCT"

# Helper function to find and save markers
find_and_save_markers <- function(seu_obj, ident_col, file_prefix) {
  Idents(seu_obj) <- ident_col
  message(glue("Finding markers for: {ident_col}"))
  
  markers <- FindAllMarkers(
    object = seu_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = FALSE
  )
  
  # Select Top 20 by LogFC
  top20 <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log2FC)
  
  write.csv(top20, file = glue("{outs_dir}/top20_genes_{file_prefix}.csv"), row.names = FALSE)
  write.csv(markers, file = glue("{outs_dir}/all_markers_{file_prefix}.csv"), row.names = FALSE)
  message(glue("Saved markers for {file_prefix}"))
}

# 3.1 Markers per Timepoint (Group)
find_and_save_markers(seu_spg, "group", "per_timepoint")

# 3.2 Markers per Stage (Early/Late)
find_and_save_markers(seu_spg, "stage", "per_stage")

# 4. DotPlots for Known Markers
# ==============================================================================
message("Generating DotPlots for known SPG markers...")

marker_genes_spg_filtered <- c("Gfra1", "Ret", "Zbtb16", "Foxo1",
                                "Kit", "Stra8", "Sycp3", "Piwil1")

# Define custom color palette
dotplot_colors <- c(rgb(252/255, 247/255, 243/255), rgb(253/255, 200/255, 180/255), 
                    rgb(229/255, 51/255, 38/255), rgb(102/255, 0/255, 13/255))

# RNA DotPlot
DefaultAssay(seu_spg) <- 'SCT'
p1 <- DotPlot(seu_spg, features = marker_genes_spg_filtered, group.by = 'celltype_spg_lvl2', 
              scale = TRUE, cluster.idents = FALSE, col.min = 0, col.max = 1.5) + 
      RotatedAxis() +
      scale_color_gradientn(colours = dotplot_colors) +
      xlab('Gene expression') + ylab('') +
      scale_y_discrete(limits = rev(levels(seu_spg$celltype_spg_lvl2))) +
      scale_x_discrete(limits = marker_genes_spg_filtered)

ggsave(filename = glue("{outs_dir}/fig_spg_supp_dotplot_rna.pdf"), plot = p1, width = 14, height = 9)

# Activity DotPlot (if assay exists)
if ("ACTIVITY" %in% names(seu_spg@assays)) {
    DefaultAssay(seu_spg) <- 'ACTIVITY'
    p2 <- DotPlot(seu_spg, features = marker_genes_spg_filtered, group.by = 'celltype_spg_lvl2', 
                  scale = TRUE, cluster.idents = FALSE, col.min = 0, col.max = 1.2) + 
          RotatedAxis() +
          scale_color_gradientn(colours = dotplot_colors) +
          xlab('Gene Activity') + ylab('') +
          scale_y_discrete(limits = rev(levels(seu_spg$celltype_spg_lvl2))) +
          scale_x_discrete(limits = marker_genes_spg_filtered)
    
    ggsave(filename = glue("{outs_dir}/fig_spg_supp_dotplot_activity.pdf"), plot = p2, width = 14, height = 9)
}

# 5. Pseudotime Visualization (RidgePlots)
# ==============================================================================
message("Visualizing Pseudotime distribution...")

# Define Timepoint Colors
timepoint_colors <- c(
  "P7" = "#1f77b4", "P14" = "#ff7f0e", "P21" = "#2ca02c",
  "P28" = "#d62728", "P35" = "#9467bd", "W10" = "#8c564b"
)

# RidgePlot by Timepoint
ridge_plot_timepoint <- RidgePlot(seu_spg, 
                                  features = "pseudotime", 
                                  group.by = "group") + 
                        scale_fill_manual(values = alpha(timepoint_colors, 0.8)) +
                        theme_minimal(base_size = 12) +
                        theme(axis.title.y = element_blank(),
                              panel.grid = element_blank())

ggsave(filename = glue("{outs_dir}/fig_spg_supp_ridgeplot_pseudotime_timepoint.pdf"), plot = ridge_plot_timepoint, width = 14, height = 9)

# 5.1 Advanced RidgePlot with Annotation Bar
# (Optional: Requires calculated bin centers and modal timepoints)
# Assuming 'annotation_data' logic is implemented if needed for complex plots.
# Here we keep the standard ridge plot which is robust.

# 6. Cellular Composition Analysis (Stacked Area Plot)
# ==============================================================================
message("Generating cellular composition stacked area plot...")

# Extract plotting data
df_for_plot <- seu_spg@meta.data[, c("pseudotime", "celltype_spg_lvl2")]

# Plot
stacked_area_plot <- ggplot(df_for_plot, aes(x = pseudotime, fill = celltype_spg_lvl2)) +
  geom_density(position = "fill", alpha = 0.8, color = "black", linewidth = 0.1) +
  theme_classic() +
  labs(title = "Cellular Composition across Pseudotime Trajectory",
       x = "Pseudotime", y = "Proportion of Cells", fill = "Cell Stage") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_brewer(palette = "Set2")

ggsave(filename = glue("{outs_dir}/spg_composition_stacked_area.pdf"), plot = stacked_area_plot, width = 10, height = 6)

message("SPG analysis visualization completed.")



marker_genes_spg_filtered <- c("Gfra1", "Ret", "Zbtb16", "Foxo1",
                                "Kit", "Stra8", "Sycp3", "Piwil1")

cell_type_col <- "celltype_spg_lvl2"
reduction_used <- "umap.rna.dims_1_25"
DefaultAssay(seu_spg) <- "SCT"

featureplot <- FeaturePlot(seu_spg, 
            features = marker_genes_spg_filtered, 
            pt.size = 0.5, 
            order = F, 
            cols = c("lightgrey", "#DE1F1F"),  # "navy",
            # cols = c("lightgrey", "navy"),  # "navy",
            # cols = c("lightgrey","#e0f2f7","#b3d8e6","#80bed9","#4aa3c8","#0f5688","#1f3368"),
            raster = F, 
            reduction = reduction_used
            # split.by = "species"
            ) & 
            xlab("") &
            ylab("") &
            theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 10), # 设置粗斜体并调整字号为16,
                   # "bold.italic"
                  legend.position = "right", 
                  # panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),  # 不要边框不过了
                  aspect.ratio = 5 / 5,
                  # 去除x轴和y轴的刻度值
                  axis.line = element_blank(),  # 隐藏轴线
                  axis.text.x = element_blank(),  # 隐藏x轴刻度标签
                  axis.text.y = element_blank(),  # 隐藏y轴刻度标签
                  axis.ticks.x = element_blank(), # 隐藏x轴刻度线
                  axis.ticks.y = element_blank()  # 隐藏y轴刻度线 
                  )

ggsave(featureplot, filename = glue("{outs_dir}/featureplot_testis_spg_filtered.pdf"), width = 20, height = 20)

