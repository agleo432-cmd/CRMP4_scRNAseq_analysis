# Load libraries
library(hdf5r)
library(Seurat)
library(dplyr)
library(tibble)
library(patchwork)
library(ggplot2)
library(Matrix)

#Set working directory (will vary by the user)
setwd("Directory")

samples <- c("GSE189070_RAW/GSM5694212_Sample_pSCI0.5d",
             "GSE189070_RAW/GSM5694213_Sample_pSCI1d",
             "GSE189070_RAW/GSM5694214_Sample_pSCI3d",
             "GSE189070_RAW/GSM5694215_Sample_pSCI7d",
             "GSE189070_RAW/GSM5694216_Sample_pSCI14d"
)

seurat_list <- list()

#
for (prefix in samples) {
  # Specify the file path
  matrix_file  <- paste0(prefix, "_matrix.mtx/matrix.mtx")
  features_file <- paste0(prefix, "_features.tsv/features.tsv")
  barcodes_file <- paste0(prefix, "_barcodes.tsv/barcodes.tsv")
  
  # Read files
  mat <- Matrix::readMM(file = matrix_file)
  features <- read.delim(features_file, header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(barcodes_file, header = FALSE, stringsAsFactors = FALSE)
  
  # Set row and column names
  gene_names <- make.unique(features$V2)
  rownames(mat) <- gene_names
  colnames(mat) <- barcodes$V1
  
  # Create a Seurat object
  seurat_obj <- CreateSeuratObject(counts = mat, project = basename(prefix),
                                   min.cells = 3, min.features = 200)
  
  # Calculate the proportion of mitochondrial genes
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-")
  
  seurat_list[[basename(prefix)]] <- seurat_obj
}

# Simplify sample names
names(seurat_list) <- basename(samples)
names(seurat_list) <- gsub("GSM[0-9]+_", "", names(seurat_list))
names(seurat_list) <- gsub("_Sample_", "_", names(seurat_list))

# Add sample name to metadata
for (i in seq_along(seurat_list)) {
  seurat_list[[i]]$sample <- names(seurat_list)[i]
}

# Merge objects (for QC visualization)
combined <- merge(seurat_list[[1]], y = seurat_list[-1])
combined$sample <- factor(combined$sample, levels = names(seurat_list))

# Violin Plot for Quality check
VlnPlot(combined,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "sample",
        ncol = 3,
        pt.size = 0.1) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  ) +
  labs(x = NULL)

# Based on the distribution observed in the violin plot obtained above,
# the thresholds for this study were defined as follows:

# 500 < nFeature_RNA < 6000
# 1000 < nCount_RNA < 40000
# percent.mt < 15