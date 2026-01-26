#!/usr/bin/env Rscript
### ==============================================================================
### 脚本1：数据预处理和质控
### 功能：数据加载、QC、过滤、标准化、降维、批次校正
### ==============================================================================

# ======== 1. 环境设置 ========
cat("\n========== 脚本1: 数据预处理和质控 ==========\n")
cat("步骤 1: 设置环境和加载包...\n")

# 设置R包路径和工作目录
.libPaths(c("/home/chengdong/R/x86_64-pc-linux-gnu-library/4.4/", .libPaths()))
setwd("/s1/chengd/project/multiome/03_brain/brain_multiome_v2_N1N2")

# 加载必要的包
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(ggplot2)
  library(patchwork)
  library(ggridges)
  library(harmony)
  library(dplyr)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(future)
})

# 设置并行处理
plan("multisession", workers = 2)
options(future.globals.maxSize = 80000 * 1024^2)

# 创建有组织的输出目录结构
cat("创建输出目录结构...\n")
dir.create("01_Preprocessing", showWarnings = FALSE)
dir.create("01_Preprocessing/QC_plots", showWarnings = FALSE)
dir.create("01_Preprocessing/Data", showWarnings = FALSE)
dir.create("01_Preprocessing/Batch_correction", showWarnings = FALSE)
dir.create("01_Preprocessing/UMAP_optimization", showWarnings = FALSE)

# ======== 2. 数据加载 ========
cat("\n步骤 2: 加载duet-seq数据...\n")

# 定义样本信息
duet_atac_samples <- c("250611DDA_N1", "250611DDA_N2")
duet_rna_samples <- gsub("DDA", "DDR", duet_atac_samples)
duet_time_points <- c("N1", "N2")

# 获取基因组注释
cat("获取基因组注释信息...\n")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

# 定义数据加载函数
load_duet_seq_data <- function(rna_path, atac_path, sample_name) {
  rna_counts <- Read10X_h5(file.path(rna_path, "filtered_feature_bc_matrix.h5"))
  atac_counts <- Read10X_h5(file.path(atac_path, "filtered_peak_bc_matrix.h5"))
  fragpath <- file.path(atac_path, "fragments.tsv.gz")
  
  # 获取RNA和ATAC的交集细胞
  common_cells <- intersect(colnames(rna_counts), colnames(atac_counts))
  cat(" -> 样本", sample_name, ": 交集细胞数 =", length(common_cells), "\n")
  
  # 创建ChromatinAssay
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts[, common_cells], 
    sep = c(":", "-"), 
    fragments = fragpath,
    genome = "mm10", 
    annotation = annotations
  )
  
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(
    counts = rna_counts[, common_cells], 
    assay = "RNA", 
    project = sample_name
  )
  seurat_obj[["ATAC"]] <- chrom_assay
  
  return(seurat_obj)
}

# 加载所有duet-seq样本
seurat_list <- list()
base_path_rna <- "/nas31/SHARE/multiome/brain/rna/force_cell_cellranger"
base_path_atac <- "/nas31/SHARE/multiome/brain/atac/force_cell_cellranger"

for(i in seq_along(duet_rna_samples)) {
  sample_rna_name <- duet_rna_samples[i]
  sample_atac_name <- duet_atac_samples[i]
  sample_path_rna <- file.path(base_path_rna, sample_rna_name, "outs")
  sample_path_atac <- file.path(base_path_atac, sample_atac_name, "outs")
  
  seurat_list[[sample_rna_name]] <- load_duet_seq_data(
    sample_path_rna, sample_path_atac, sample_rna_name
  )
  
  # 添加元数据
  seurat_list[[sample_rna_name]]$orig.ident <- sample_rna_name
  seurat_list[[sample_rna_name]]$platform <- "duet_seq"
  seurat_list[[sample_rna_name]]$replicate <- paste0("rep_", i)
  seurat_list[[sample_rna_name]]$time_point <- duet_time_points[i]
  
  cat("成功加载:", sample_rna_name, "\n")
}

# 合并数据
cat("\n合并所有样本...\n")
if(length(seurat_list) > 1) {
  merged_seurat <- merge(seurat_list[[1]], seurat_list[-1], 
                         add.cell.ids = names(seurat_list))
} else {
  merged_seurat <- seurat_list[[1]]
}

total_cells_initial <- ncol(merged_seurat)
cat("✓ 合并后总细胞数:", total_cells_initial, "\n")
cat("✓ 样本分布:\n")
print(table(merged_seurat$orig.ident))

# ======== 3. 质控计算 ========
cat("\n步骤 3: 计算质控指标...\n")

# RNA QC指标
DefaultAssay(merged_seurat) <- "RNA"
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^mt-")

# ATAC QC指标
DefaultAssay(merged_seurat) <- "ATAC"
merged_seurat <- NucleosomeSignal(merged_seurat)
plan("sequential")  # TSS计算需要单线程
merged_seurat <- TSSEnrichment(merged_seurat, fast = TRUE)
plan("multisession", workers = 2)

# 计算FRiP
cat("计算FRiP指标...\n")
calculate_frip <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "ATAC"
  frag_list <- Fragments(seurat_obj[['ATAC']])
  
  if (length(frag_list) == 0) {
    stop("错误: 未找到fragment文件")
  }
  
  seurat_obj$total_fragments <- NA
  
  for (i in seq_along(frag_list)) {
    frag_path <- frag_list[[i]]@path
    cat("  处理fragment文件", i, ":", basename(frag_path), "\n")
    
    if (!file.exists(frag_path)) {
      warning("文件不存在: ", frag_path)
      next
    }
    
    tryCatch({
      cells_in_frag <- frag_list[[i]]@cells
      if(length(cells_in_frag) == 0) next
      
      seurat_cells <- names(cells_in_frag)
      original_barcodes <- unname(cells_in_frag)
      
      file_counts <- CountFragments(
        fragments = frag_path,
        cells = original_barcodes,
        verbose = FALSE
      )
      
      counts_lookup <- setNames(file_counts$frequency_count, file_counts$CB)
      
      for(j in seq_along(seurat_cells)) {
        seurat_cell_name <- seurat_cells[j]
        original_barcode <- cells_in_frag[seurat_cell_name]
        
        if(original_barcode %in% names(counts_lookup)) {
          seurat_obj$total_fragments[seurat_cell_name] <- counts_lookup[original_barcode]
        }
      }
    }, error = function(e) {
      warning("处理文件时出错: ", e$message)
    })
  }
  
  # 对未计算的细胞进行估算
  unmapped_cells <- is.na(seurat_obj$total_fragments)
  if (any(unmapped_cells)) {
    seurat_obj$total_fragments[unmapped_cells] <- seurat_obj$nCount_ATAC[unmapped_cells] * 2
  }
  
  # 计算FRiP
  peak_fragments <- seurat_obj$nCount_ATAC / 2
  frip_values <- peak_fragments / seurat_obj$total_fragments
  frip_values[is.infinite(frip_values) | is.na(frip_values)] <- 0
  frip_values[frip_values > 1] <- 1
  seurat_obj$FRiP <- frip_values
  
  return(seurat_obj)
}

merged_seurat <- calculate_frip(merged_seurat)

# ======== 4. QC可视化（过滤前）========
cat("\n步骤 4: 生成过滤前QC图...\n")

# RNA QC图
p1 <- VlnPlot(merged_seurat, features = "nFeature_RNA", group.by = "orig.ident", pt.size = 0) + 
  scale_y_log10() + NoLegend() + theme(axis.title.x = element_blank())
p2 <- VlnPlot(merged_seurat, features = "nCount_RNA", group.by = "orig.ident", pt.size = 0) + 
  scale_y_log10() + NoLegend() + theme(axis.title.x = element_blank())
p3 <- VlnPlot(merged_seurat, features = "percent.mt", group.by = "orig.ident", pt.size = 0) + 
  NoLegend() + theme(axis.title.x = element_blank())
rna_qc_pre <- p1 | p2 | p3
rna_qc_pre <- rna_qc_pre + plot_annotation(title = "RNA QC Metrics (Pre-filter)")
ggsave("01_Preprocessing/QC_plots/RNA_QC_pre_filter.pdf", rna_qc_pre, width = 15, height = 5)

# ATAC QC图
p4 <- VlnPlot(merged_seurat, features = "nFeature_ATAC", group.by = "orig.ident", pt.size = 0) + 
  scale_y_log10() + NoLegend() + theme(axis.title.x = element_blank())
p5 <- VlnPlot(merged_seurat, features = "nCount_ATAC", group.by = "orig.ident", pt.size = 0) + 
  scale_y_log10() + NoLegend() + theme(axis.title.x = element_blank())
p6 <- VlnPlot(merged_seurat, features = "nucleosome_signal", group.by = "orig.ident", pt.size = 0) + 
  NoLegend() + theme(axis.title.x = element_blank())
p7 <- VlnPlot(merged_seurat, features = "TSS.enrichment", group.by = "orig.ident", pt.size = 0) + 
  NoLegend() + theme(axis.title.x = element_blank())
p8 <- VlnPlot(merged_seurat, features = "FRiP", group.by = "orig.ident", pt.size = 0) + 
  NoLegend() + theme(axis.title.x = element_blank())
atac_qc_pre <- p4 | p5 | p6 | p7 | p8
atac_qc_pre <- atac_qc_pre + plot_annotation(title = "ATAC QC Metrics (Pre-filter)")
ggsave("01_Preprocessing/QC_plots/ATAC_QC_pre_filter.pdf", atac_qc_pre, width = 25, height = 5)

# ======== 5. 细胞过滤 ========
cat("\n步骤 5: 执行细胞过滤...\n")
cat("过滤阈值:\n")
cat("  nFeature_RNA: 500-9000\n")
cat("  percent.mt: < 20\n")
cat("  nCount_ATAC: > 1500\n")
cat("  TSS.enrichment: > 2\n")
cat("  FRiP: > 0.20\n")
cat("  nucleosome_signal: < 1.5\n")

merged_seurat <- subset(
  x = merged_seurat,
  subset = nFeature_RNA > 500 &
    nFeature_RNA < 9000 &
    percent.mt < 20 &
    nCount_ATAC > 1500 &
    TSS.enrichment > 2 &
    FRiP > 0.20 &
    nucleosome_signal < 1.5
)

cells_after_qc <- ncol(merged_seurat)
cat("QC过滤后细胞数:", cells_after_qc, 
    "(保留率:", round(cells_after_qc/total_cells_initial*100, 1), "%)\n")

# QC可视化（过滤后）
cat("生成过滤后QC图...\n")

# RNA QC图（过滤后）
p1_post <- VlnPlot(merged_seurat, features = "nFeature_RNA", group.by = "orig.ident", pt.size = 0) + 
  scale_y_log10() + NoLegend() + theme(axis.title.x = element_blank())
p2_post <- VlnPlot(merged_seurat, features = "nCount_RNA", group.by = "orig.ident", pt.size = 0) + 
  scale_y_log10() + NoLegend() + theme(axis.title.x = element_blank())
p3_post <- VlnPlot(merged_seurat, features = "percent.mt", group.by = "orig.ident", pt.size = 0) + 
  NoLegend() + theme(axis.title.x = element_blank())
rna_qc_post <- p1_post | p2_post | p3_post
rna_qc_post <- rna_qc_post + plot_annotation(title = "RNA QC Metrics (Post-filter)")
ggsave("01_Preprocessing/QC_plots/RNA_QC_post_filter.pdf", rna_qc_post, width = 15, height = 5)

# ATAC QC图（过滤后）
p4_post <- VlnPlot(merged_seurat, features = "nFeature_ATAC", group.by = "orig.ident", pt.size = 0) + 
  scale_y_log10() + NoLegend() + theme(axis.title.x = element_blank())
p5_post <- VlnPlot(merged_seurat, features = "nCount_ATAC", group.by = "orig.ident", pt.size = 0) + 
  scale_y_log10() + NoLegend() + theme(axis.title.x = element_blank())
p6_post <- VlnPlot(merged_seurat, features = "nucleosome_signal", group.by = "orig.ident", pt.size = 0) + 
  NoLegend() + theme(axis.title.x = element_blank())
p7_post <- VlnPlot(merged_seurat, features = "TSS.enrichment", group.by = "orig.ident", pt.size = 0) + 
  NoLegend() + theme(axis.title.x = element_blank())
p8_post <- VlnPlot(merged_seurat, features = "FRiP", group.by = "orig.ident", pt.size = 0) + 
  NoLegend() + theme(axis.title.x = element_blank())
atac_qc_post <- p4_post | p5_post | p6_post | p7_post | p8_post
atac_qc_post <- atac_qc_post + plot_annotation(title = "ATAC QC Metrics (Post-filter)")
ggsave("01_Preprocessing/QC_plots/ATAC_QC_post_filter.pdf", atac_qc_post, width = 25, height = 5)

# ======== 6. 双细胞检测 ========
cat("\n步骤 6: 双细胞检测与过滤...\n")
seurat_list_dbl <- SplitObject(merged_seurat, split.by = "orig.ident")

plan("sequential")
seurat_list_dbl <- lapply(names(seurat_list_dbl), function(sample_name) {
  obj <- seurat_list_dbl[[sample_name]]
  cat("  处理样本:", sample_name, "\n")
  
  tryCatch({
    sce <- as.SingleCellExperiment(obj, assay = "RNA")
    sce <- scDblFinder(sce, verbose = FALSE)
    obj$scDblFinder.score <- sce$scDblFinder.score
    obj$scDblFinder.class <- sce$scDblFinder.class
    return(obj)
  }, error = function(e) {
    warning("双细胞检测失败: ", e$message)
    obj$scDblFinder.score <- 0
    obj$scDblFinder.class <- "singlet"
    return(obj)
  })
})
names(seurat_list_dbl) <- names(SplitObject(merged_seurat, split.by = "orig.ident"))
plan("multisession", workers = 2)

merged_seurat <- merge(seurat_list_dbl[[1]], y = seurat_list_dbl[-1])

# 双细胞检测可视化
p_doublet_score <- VlnPlot(merged_seurat, features = "scDblFinder.score", 
                           group.by = "orig.ident", pt.size = 0) +
  labs(title = "Doublet Scores by Sample", y = "Doublet Score")

p_doublet_class <- merged_seurat@meta.data %>%
  ggplot(aes(x = orig.ident, fill = scDblFinder.class)) +
  geom_bar(position = "fill") +
  labs(title = "Doublet Classification by Sample", 
       x = "Sample", y = "Proportion", fill = "Classification") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

doublet_combined <- p_doublet_score / p_doublet_class
ggsave("01_Preprocessing/QC_plots/Doublet_Detection.pdf", doublet_combined, width = 12, height = 8)

# 过滤双细胞
total_doublets <- sum(merged_seurat$scDblFinder.class == "doublet")
cat("检测到双细胞数:", total_doublets, "\n")
merged_seurat <- subset(merged_seurat, subset = scDblFinder.class == "singlet")
cells_after_doublet <- ncol(merged_seurat)
cat("双细胞过滤后细胞数:", cells_after_doublet, "\n")

# ======== 7. 标准化与降维 ========
cat("\n步骤 7: 标准化与降维...\n")

# RNA标准化
DefaultAssay(merged_seurat) <- "RNA"
merged_seurat <- SCTransform(merged_seurat, vars.to.regress = "percent.mt", verbose = FALSE)
merged_seurat <- RunPCA(merged_seurat, verbose = FALSE, npcs = 50)

# ATAC标准化
DefaultAssay(merged_seurat) <- "ATAC"
merged_seurat <- RunTFIDF(merged_seurat, method = 1)
merged_seurat <- FindTopFeatures(merged_seurat, min.cutoff = "q5")
merged_seurat <- RunSVD(merged_seurat, n = 50, verbose = FALSE)

# Elbow plots
p_elbow_pca <- ElbowPlot(merged_seurat, ndims = 50, reduction = "pca") + 
  ggtitle("Elbow Plot - RNA PCA")
ggsave("01_Preprocessing/QC_plots/Elbow_Plot_PCA.pdf", p_elbow_pca, width = 8, height = 6)

p_elbow_lsi <- ElbowPlot(merged_seurat, ndims = 50, reduction = "lsi") + 
  ggtitle("Elbow Plot - ATAC LSI")
ggsave("01_Preprocessing/QC_plots/Elbow_Plot_LSI.pdf", p_elbow_lsi, width = 8, height = 6)

# ======== 8. Harmony批次校正 ========
cat("\n步骤 8: Harmony批次校正...\n")

# RNA Harmony
DefaultAssay(merged_seurat) <- "SCT"
merged_seurat <- RunHarmony(
  object = merged_seurat, 
  group.by.vars = c("orig.ident", "replicate"),
  reduction = "pca", 
  assay.use = "SCT", 
  reduction.save = "harmony_rna",
  verbose = FALSE
)

# ATAC Harmony
DefaultAssay(merged_seurat) <- "ATAC"
merged_seurat <- RunHarmony(
  object = merged_seurat, 
  group.by.vars = c("orig.ident", "replicate"),
  reduction = "lsi", 
  assay.use = "ATAC", 
  project.dim = FALSE, 
  reduction.save = "harmony_atac",
  dims.use = 2:50,
  verbose = FALSE
)

# 批次校正效果可视化
pca_before <- DimPlot(merged_seurat, reduction = "pca", group.by = "orig.ident") + 
  ggtitle("RNA PCA before Harmony")
pca_after <- DimPlot(merged_seurat, reduction = "harmony_rna", group.by = "orig.ident") + 
  ggtitle("RNA PCA after Harmony")
pca_comparison <- pca_before | pca_after
ggsave("01_Preprocessing/Batch_correction/RNA_Harmony_comparison.pdf", pca_comparison, width = 16, height = 7)

lsi_before <- DimPlot(merged_seurat, reduction = "lsi", group.by = "orig.ident", dims = c(2,3)) + 
  ggtitle("ATAC LSI before Harmony")
lsi_after <- DimPlot(merged_seurat, reduction = "harmony_atac", group.by = "orig.ident", dims = c(2,3)) + 
  ggtitle("ATAC LSI after Harmony")
lsi_comparison <- lsi_before | lsi_after
ggsave("01_Preprocessing/Batch_correction/ATAC_Harmony_comparison.pdf", lsi_comparison, width = 16, height = 7)

# ======== 9. UMAP参数优化 ========
cat("\n步骤 9: UMAP参数优化...\n")

# RNA UMAP参数测试
rna_dims_options <- list(
  "dims_1_20" = 1:20, "dims_1_25" = 1:25, "dims_1_30" = 1:30,
  "dims_1_35" = 1:35, "dims_1_40" = 1:40, "dims_1_50" = 1:50
)

for(dim_name in names(rna_dims_options)) {
  dims <- rna_dims_options[[dim_name]]
  merged_seurat <- RunUMAP(
    merged_seurat,
    reduction = "harmony_rna",
    dims = dims,
    reduction.name = paste0("rna.umap.", dim_name),
    reduction.key = paste0("rnaUMAP", gsub("dims_", "", dim_name), "_"),
    n.neighbors = 30,
    min.dist = 0.3,
    spread = 1,
    verbose = FALSE
  )
}

# ATAC UMAP参数测试
atac_dims_options <- list(
  "dims_2_20" = 2:20, "dims_2_25" = 2:25, "dims_2_30" = 2:30,
  "dims_2_35" = 2:35, "dims_2_40" = 2:40, "dims_2_49" = 2:49
)

for(dim_name in names(atac_dims_options)) {
  dims <- atac_dims_options[[dim_name]]
  merged_seurat <- RunUMAP(
    merged_seurat, 
    reduction = "harmony_atac", 
    dims = dims, 
    reduction.name = paste0("atac.umap.", dim_name), 
    reduction.key = paste0("atacUMAP", gsub("dims_", "", dim_name), "_"),
    n.neighbors = 30,
    min.dist = 0.4,
    spread = 1,
    verbose = FALSE
  )
}

# 可视化UMAP参数对比
rna_plots <- lapply(names(rna_dims_options), function(dim_name) {
  DimPlot(merged_seurat, reduction = paste0("rna.umap.", dim_name), 
          group.by = "orig.ident", shuffle = TRUE) +
    ggtitle(paste("RNA UMAP", gsub("dims_", "dims ", dim_name))) +
    NoLegend() + theme_minimal()
})
combined_rna_umaps <- wrap_plots(rna_plots, ncol = 3)
ggsave("01_Preprocessing/UMAP_optimization/RNA_UMAP_dimensions.pdf", 
       combined_rna_umaps, width = 18, height = 12)

atac_plots <- lapply(names(atac_dims_options), function(dim_name) {
  DimPlot(merged_seurat, reduction = paste0("atac.umap.", dim_name), 
          group.by = "orig.ident", shuffle = TRUE) +
    ggtitle(paste("ATAC UMAP", gsub("dims_", "dims ", dim_name))) +
    NoLegend() + theme_minimal()
})
combined_atac_umaps <- wrap_plots(atac_plots, ncol = 3)
ggsave("01_Preprocessing/UMAP_optimization/ATAC_UMAP_dimensions.pdf", 
       combined_atac_umaps, width = 18, height = 12)

# ======== 10. 保存结果 ========
cat("\n步骤 10: 保存预处理结果...\n")
saveRDS(merged_seurat, file = "01_Preprocessing/Data/preprocessed_seurat.rds")

# 生成预处理总结
preprocessing_summary <- data.frame(
  Step = c("Initial cells", "After QC filter", "After doublet removal", "Final cells"),
  Cell_count = c(total_cells_initial, cells_after_qc, cells_after_doublet, ncol(merged_seurat)),
  Retention_rate = c(100, 
                     round(cells_after_qc/total_cells_initial*100, 1),
                     round(cells_after_doublet/total_cells_initial*100, 1),
                     round(ncol(merged_seurat)/total_cells_initial*100, 1))
)

write.csv(preprocessing_summary, "01_Preprocessing/preprocessing_summary.csv", row.names = FALSE)

cat("\n========== 预处理完成总结 ==========\n")
print(preprocessing_summary)
cat("\n输出文件已保存至 01_Preprocessing/ 文件夹\n")
cat("下一步: 运行脚本2进行聚类和细胞类型注释\n")
cat("=====================================\n")
