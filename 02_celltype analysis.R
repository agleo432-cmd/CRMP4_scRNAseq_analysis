# This script inputs 6 files (post SCI 0.5d, 1d, 3d, 7d, 14d)
# It preforms QC, combines and normalizes the files, then assigns
# identity to clusters based on marker genes.
# Finally, we searched the Dpysl3 expression in all celltpes using featureplot

# Load libraries
library(hdf5r)
library(Seurat)
library(dplyr)
library(tibble)
library(patchwork)
library(ggplot2)
library(Matrix)
library(clustree)

#set working directory (will vary by the user)
setwd("directory")

#########   Input Data    ###########################################

seurat_list <- list()

samples <- c("GSE189070_RAW/GSM5694212_Sample_pSCI0.5d",
             "GSE189070_RAW/GSM5694213_Sample_pSCI1d",
             "GSE189070_RAW/GSM5694214_Sample_pSCI3d",
             "GSE189070_RAW/GSM5694215_Sample_pSCI7d",
             "GSE189070_RAW/GSM5694216_Sample_pSCI14d"
)

#Read file and QC
for (prefix in samples) {
  
  message("Processing: ", basename(prefix))
  
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
  
  seurat_obj <- RenameCells(seurat_obj, add.cell.id = basename(prefix))
  
# Normalize
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 5000)
  
  seurat_obj$sample <- basename(prefix)
  
  seurat_list[[basename(prefix)]] <- seurat_obj
}

# Cell count for each sample (initial)
sapply(seurat_list, ncol)

# Cell count will be as follows:
#GSM5694212_Sample_pSCI0.5d   GSM5694213_Sample_pSCI1d 
#15951                      11309 
#GSM5694214_Sample_pSCI3d   GSM5694215_Sample_pSCI7d 
#11729                      15808 
#GSM5694216_Sample_pSCI14d 
#9301 

# Number of genes per sample
sapply(seurat_list, nrow)

# Number of genes will be as follows:
#GSM5694212_Sample_pSCI0.5d   GSM5694213_Sample_pSCI1d 
#17299                      17775 
#GSM5694214_Sample_pSCI3d   GSM5694215_Sample_pSCI7d 
#16940                      17480 
#GSM5694216_Sample_pSCI14d 
#16938 

#set seed for ensuring reproducibility
set.seed(123)

# To reduce data volume, set the number of cells per sample to 5000.
max_cells <- 5000  
for (i in 1:length(seurat_list)) {
  seurat_obj <- seurat_list[[i]]
  if (ncol(seurat_obj) > max_cells) {
    sampled_cells <- sample(colnames(seurat_obj), max_cells)
    seurat_list[[i]] <- subset(seurat_obj, cells = sampled_cells)
  }
}

# Cell count for each sample (after)
sapply(seurat_list, ncol)

# Number of genes will be as follows:
#GSM5694212_Sample_pSCI0.5d   GSM5694213_Sample_pSCI1d 
#5000                      5000 
#GSM5694214_Sample_pSCI3d   GSM5694215_Sample_pSCI7d 
#5000                      5000 
#GSM5694216_Sample_pSCI14d 
#5000

#Finding anchors for integration
min_cells <- min(sapply(seurat_list, ncol))

anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  dims = 1:30,
  k.filter = min_cells
)

#Merge sample data
integrated_data <- IntegrateData(anchorset = anchors, dims = 1:30)

#Scaling、PCA、UMAP
DefaultAssay(integrated_data) <- "integrated"  

integrated_data <- ScaleData(integrated_data, verbose = FALSE)
integrated_data <- RunPCA(integrated_data, verbose = FALSE)
integrated_data <- RunUMAP(integrated_data, dims = 1:30)

# Confirm UMAP (classified by "Day per injury")
DimPlot(integrated_data, group.by = "sample") +
  coord_fixed() +
  theme(
    plot.title = element_blank()
  )

# Save integrated data as rds
saveRDS(integrated_data, file = "integrated_data_GSE189070.rds")

#########   Clustering and Annotation    ##########################

integrated_data <- FindNeighbors(integrated_data,dims=1:30)
integrated_data <- FindClusters(integrated_data, resolution = c(0.2, 0.4, 0.6, 1.0))

# Make clustree
clustree(integrated_data@meta.data, prefix = "integrated_snn_res.")

# Set resolution as 0.6
Idents(integrated_data) <- "integrated_snn_res.0.6"

current_ids <- levels(Idents(integrated_data))

new_ids <- sort(as.numeric(current_ids))

Idents(integrated_data) <- factor(
  Idents(integrated_data),
  levels = as.character(new_ids)
)

table(Idents(integrated_data))
# 19 clusters identified

# Visualize clusters on UMAP
DimPlot(integrated_data, reduction = "umap", label = TRUE, pt.size = 0.5) +
  coord_fixed() +
  theme(
    plot.title = element_blank()
  )

#Strategy for clustering
# 1. Perform coarse clustering using feature plots of representative 
#    markers for astrocytes, microglia, and ODC.
# 2. Extract the top markers from each cluster and identify cell types
#    from marker types.

FeaturePlot(
  integrated_data,
  features = c("Aif1", "Fcrls", "P2ry12", "Cx3cr1", "Gfap", "Aldh1l1", "Aqp4", "Slc1a3", "S100b","Cspg4", "Olig1", "Mbp"),
  ncol = 4
) +
  coord_fixed()

# Marker Gene Search for Each Cluster
# Extract genes with logFC ≥ 0.25 and expressed in ≥ 25% of the cluster
markers <- FindAllMarkers(
  integrated_data,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Extract the top 10 (5) markers from each cluster
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top5 <- markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

# Save marker data (Uploaded to Github)
write.csv(markers, "markers_res0.6_GSE189070.csv", row.names = FALSE)

# Confirm with heatmap
DoHeatmap(integrated_data, features = top5$gene) +
  scale_fill_gradientn(colors = c("navy", "white", "firebrick3"))

# We annotated the clusters by referring both the csv file and heatmap
# comparing with the clustered UMAP (resolution: 0.6).
# e.g. From the feature plot, cluster 0 has a pattern of expression 
# similar to microglia. From the table of the top 10 markers, Ch25h
# and Ccl3, which have been reported to promote neuroinflammatory
# responses, were detected, so we defined it as M1 microglia.

Idents(integrated_data)

#Labeling each cluster

cluster_to_celltype <- c(
  "0" = "M1 microglia",
  "1" = "M2 microglia",
  "2" = "Microglia",
  "3" = "Microglia",
  "4" = "Granulocyte",
  "5" = "Microglia",
  "6" = "ODC",
  "7" = "Microglia",
  "8" = "Macrophage",
  "9" = "Microglia",
  "10" = "Dendritic cell",
  "11" = "B cell",
  "12" = "Neutrophil",
  "13" = "OPC",
  "14" = "NK cells",
  "15" = "Astrocyte",
  "16" = "Astrocyte",
  "17" = "Red blood cell",
  "18" = "Neutrophil precursor cells",
  "19" = "pDC"
)

integrated_data <- RenameIdents(
  integrated_data,
  cluster_to_celltype
)

integrated_data$celltype <- Idents(integrated_data)

table(integrated_data$celltype)

# gain clustering results
cluster_counts <- table(Idents(integrated_data))

total_cells <- sum(cluster_counts)

relative_frequency <- cluster_counts / total_cells * 100

cluster_data <- data.frame(
  Cluster = names(cluster_counts),
  CellCount = as.numeric(cluster_counts),
  RelativeFrequency = relative_frequency
)

# Cluster sorting according to relative frequency
cluster_data$Cluster <- factor(cluster_data$Cluster, 
                               levels = cluster_data$Cluster[order(cluster_data$RelativeFrequency.Freq, decreasing = FALSE)])

ggplot(cluster_data, aes(x = Cluster, y = RelativeFrequency.Freq, fill = RelativeFrequency.Freq)) +
  geom_bar(stat = "identity", width = 0.7) +  
  coord_flip() +  
  scale_fill_gradient(low = "lightblue", high = "darkblue") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, size = 12),  
    axis.text.y = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold"),  
    plot.margin = margin(10, 30, 10, 10),  
    legend.position = c(0.85, 0.1),  
    legend.title = element_text(size = 12),  
    legend.text = element_text(size = 10) 
  ) +  
  labs(y = "Relative Frequency (%)", title = "Cluster Relative Frequency")+
  geom_text(aes(label = paste0(CellCount, " (", round(RelativeFrequency.Freq, 2), "%)")), 
            position = position_nudge(x = 0.4), 
            color = "black", size = 4, hjust = 0)

#########   Identification of CRMP4 expressing cells    ###############

FeaturePlot(
  integrated_data,
  features = c("Dpysl3"),
  min.cutoff = "q2",  
  max.cutoff = "q100",  
  cols = c("gray", "blue"),  
  pt.size = 1.0  
) + 
  coord_fixed()