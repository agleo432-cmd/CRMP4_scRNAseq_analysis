# Single-cell RNA-seq data (GSE189070) were processed in 
# Seurat, and astrocytes were selected based on expression
# of astrocyte markers and absence of immune markers.
# Filtered cells from each time point were normalized and integrated using 
# anchor-based integration.
# Clustering, UMAP visualization, cell cycle scoring, and 
# Dpysl3 (CRMP4) expression analyses were then performed to characterize 
# astrocyte subpopulations.

library(hdf5r)
library(Seurat)
library(dplyr)
library(tibble)
library(patchwork)
library(ggplot2)
library(Matrix)
library(clusterProfiler)
library(org.Mm.eg.db)
library(clustree)
library(ggrepel)

#Set working directory (will vary by the user)
setwd("directory")

#########   Input Data    ###########################################

samples <- c("GSE189070_RAW/GSM5694212_Sample_pSCI0.5d",
             "GSE189070_RAW/GSM5694213_Sample_pSCI1d",
             "GSE189070_RAW/GSM5694214_Sample_pSCI3d",
             "GSE189070_RAW/GSM5694215_Sample_pSCI7d",
             "GSE189070_RAW/GSM5694216_Sample_pSCI14d"
)

# marker genes for screening astrocytes
astro_markers <- c("Gfap", "Gss", "Aldh1l1", "Slc1a3", "AldoC", "Aqp4", "Id3", "Fabp7", "Lcn2")
immune_markers <- c("Dock2", "Cx3cr1")

seurat_list <- list()

#Read file and QC
for (prefix in samples) {
  
  matrix_file  <- file.path(paste0(prefix, "_matrix.mtx"), "matrix.mtx")
  features_file <- file.path(paste0(prefix, "_features.tsv"), "features.tsv")
  barcodes_file <- file.path(paste0(prefix, "_barcodes.tsv"), "barcodes.tsv")
  
  mat <- Matrix::readMM(file = matrix_file)
  features <- read.delim(features_file, header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(barcodes_file, header = FALSE, stringsAsFactors = FALSE)
  
  gene_names <- make.unique(features$V2)
  rownames(mat) <- gene_names
  colnames(mat) <- barcodes$V1
  
  seurat_obj <- CreateSeuratObject(counts = mat, project = basename(prefix),
                                   min.cells = 3, min.features = 200)
  
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  # Setting thresholds set by QC analysis
  seurat_obj <- subset(seurat_obj,
                       subset = nFeature_RNA > 500 &
                         nFeature_RNA < 6000 &
                         nCount_RNA > 1000 &
                         nCount_RNA < 40000 &
                         percent.mt < 15)
  
  # Normalize
  seurat_obj <- NormalizeData(seurat_obj)
  
  expr_mat <- GetAssayData(seurat_obj, slot = "data")
  
  # Condition 1: One of the Astrocyte markers is expressing
  astro_any <- Matrix::colSums(expr_mat[intersect(astro_markers, rownames(expr_mat)), , drop=FALSE] > 0) > 0
  
  # Condition 2: Immune cell relative markers are not expressing
  has_Dock2 <- "Dock2" %in% rownames(expr_mat)
  has_Cx3cr1 <- "Cx3cr1" %in% rownames(expr_mat)
  
  immu_neg <- rep(TRUE, ncol(expr_mat))
  if (has_Dock2) immu_neg <- immu_neg & (expr_mat["Dock2", ] == 0)
  if (has_Cx3cr1) immu_neg <- immu_neg & (expr_mat["Cx3cr1", ] == 0)
  
  # extract cells which meet both condition
  is_astro <- astro_any & immu_neg
  
  seu_astro <- seurat_obj[, is_astro]
  
  message(basename(prefix), " astro cells: ", ncol(seu_astro))
  
  seurat_list[[basename(prefix)]] <- seu_astro
}

# integrate data
for (i in 1:length(seurat_list)) {
  seurat_list[[i]] <- NormalizeData(seurat_list[[i]])
  seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]])
}

anchors <- FindIntegrationAnchors(object.list = seurat_list)
astro_merged <- IntegrateData(anchorset = anchors)

# save integrated seurat data (rds uploaded to Github)
saveRDS(astro_merged, file = "integrated_data_astro_injured_GSE189070.rds")

#########   Clustering    ###########################################

DefaultAssay(astro_merged) <- "integrated"

astro_merged <- ScaleData(astro_merged)
astro_merged <- RunPCA(astro_merged)
astro_merged <- RunUMAP(astro_merged, dims = 1:20)

astro_merged <- FindNeighbors(astro_merged,dims=1:30)
astro_merged <- FindClusters(astro_merged, resolution = c(0.2, 0.4, 0.6, 1.0))

clustree(astro_merged@meta.data, prefix = "integrated_snn_res.")


#feature plots for identifying the feature of each cluster (+ CRMP4(Dpysl3))
features <- c("Dpysl3","Gfap","Aldh1l1","Sox10","Olig1","Cryab", "Ccnb1", "Cdk1", "Mki67")

plots <- FeaturePlot(
  astro_merged,
  features = features,
  min.cutoff = "q2",
  max.cutoff = "q100",
  cols = c("lightgray", "blue"),
  combine = FALSE
)

plots <- lapply(plots, function(p) {
  p + 
    coord_fixed() +
    theme(
      aspect.ratio = 1,
      plot.margin = margin(2, 2, 2, 2)
    )
})

#Fig 5b
wrap_plots(plots, ncol = 3)

Idents(astro_merged) <- "integrated_snn_res.0.4"

#Fig 5a
DimPlot(astro_merged, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP (Resolution 0.4)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14)
  ) +
  coord_fixed()


# cell cycle phase marker sets
s.genes <- c("Mcm5", "Tyms", "Mcm2", "Mcm4", "Rrm1", "Gins2",
             "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1", "Hells",
             "Nasp", "Rad51ap1", "Gmnn", "Slbp", "Ccne2",
             "Atad2", "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin",
             "Clspn")

g2m.genes <- c("Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", "Ndc80",
               "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo", "Cenpf", "Tacc3", "Fam64a",
               "Smc4", "Ccnb2", "Ckap2l", "Ckap2", "Aurkb", "Kif11",
               "Tubb4b", "Kif20b", "Cdca3", "Cdc20", "Dlgap5", "Cdca2", "Cdca8", "Ect2", "Kif23",
               "Hmmr", "Aurka", "Anln", "Cenpe", "Cenpa")

#cell cycle phase scoring
astro_merged <- CellCycleScoring(astro_merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

p1 <- DimPlot(astro_merged, group.by = "Phase", label = FALSE) +
  ggtitle("Cell Cycle Phase") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.text  = element_text(size = 14)
  )

#Fig 5c
p1

# Comfirm CRMP4 expression in each cluster (Fig 5d)
VlnPlot(astro_merged, features = "Dpysl3") +
  theme(
    plot.title   = element_text(size = 18, hjust = 0.5),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  )