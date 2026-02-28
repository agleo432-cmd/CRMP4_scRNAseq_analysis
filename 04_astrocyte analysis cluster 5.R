# This script identifies CRMP4 (Dpysl3)–positive and –negative cells within 
# astrocyte clusters 5, performs differential gene expression 
# analysis between these two states, and visualizes the results using a 
# volcano plot. Gene Ontology (BP, CC, MF) and KEGG pathway enrichment analyses
# are then conducted to characterize functional differences associated with 
# Dpysl3 expression.

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

wrap_plots(plots, ncol = 3)

Idents(astro_merged) <- "integrated_snn_res.0.4"

DimPlot(astro_merged, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP (Resolution 0.4)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14)
  ) +
  coord_fixed()

##### Expression & function analysis of CRMP4 expressing astrocytes (Cluster 5)  #####

# Define CRMP4 (Dpysl3) expression status in subset clusters

# Comfirm CRMP4 expression in each cluster
VlnPlot(astro_merged, features = "Dpysl3") +
  theme(
    plot.title   = element_text(size = 18, hjust = 0.5),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  )

# extract cluster 5 astrocytes
cluster5 <- subset(astro_merged, idents = "5")

gene_name <- "Dpysl3"
expr_mat <- GetAssayData(cluster5, slot = "data")

Dpysl3_pos <- expr_mat[gene_name, ] > 0

cluster5$Dpysl3_status <- ifelse(Dpysl3_pos, "Dpysl3_pos", "Dpysl3_neg")

table(cluster5$Dpysl3_status)

# Differential gene expression analysis (Cluster 5)

deg_cl5 <- FindMarkers(
  astro_merged,
  ident.1 = 5,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

head(deg_cl5)

# save as CSV
write.csv(deg_cl5, file = "DEG cluster5 hemi SCI astrocyte.csv")

deg_cl5$log10_pval_plot <- pmin(
  -log10(deg_cl5$p_val_adj),
  300
)

# Annotate DEG directionality
deg_cl5$group <- "non-significant"
deg_cl5$group[deg_cl5$avg_log2FC > 0.25 & deg_cl5$p_val_adj < 0.05] <- "Up (Cluster5)"
deg_cl5$group[deg_cl5$avg_log2FC < -0.25 & deg_cl5$p_val_adj < 0.05] <- "Down"

deg_cl5$gene <- rownames(deg_cl5)

ggplot(deg_cl5, aes(x = avg_log2FC, y = log10_pval_plot, color = group)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c(
    "Up (Cluster5)" = "red",
    "Down" = "blue",
    "non-significant" = "gray"
  )) +
  theme_bw(base_size = 14) +
  labs(
    title = "Volcano Plot: Cluster 5",
    x = "log2 Fold Change",
    y = "-log10(Adjusted P-value)",
    color = "Group"
  ) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = subset(deg_cl5, group != "non-significant"),
    aes(label = gene),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3
  )




# Differential gene expression analysis (Dpysl3+ vs Dpysl3−)

Idents(cluster5) <- "Dpysl3_status"

# Detecting DEGs
Dpysl3_DEG <- FindMarkers(cluster5,
                          ident.1 = "Dpysl3_pos",
                          ident.2 = "Dpysl3_neg",
                          logfc.threshold = 0.25,
                          min.pct = 0.1,
                          test.use = "wilcox") 

Dpysl3_DEG %>%
  rownames_to_column("gene") %>%
  arrange(p_val_adj) %>%
  head(20)


# Volcano plot visualization

deg <- Dpysl3_DEG

# save as CSV
write.csv(deg, file = "DEG_cluster5 hemi SCI astrocyte CRMP4+.csv")

deg$log10_pval <- -log10(deg$p_val_adj + 1e-300)

# Annotate DEG directionality
deg$group <- "non-significant"
deg$group[deg$avg_log2FC > 0.25 & deg$p_val_adj < 0.05] <- "Up (Dpysl3⁺)"
deg$group[deg$avg_log2FC < -0.25 & deg$p_val_adj < 0.05] <- "Down (Dpysl3⁻)"

deg$gene <- rownames(deg)

# Generate volcano plot
ggplot(deg, aes(x = avg_log2FC, y = log10_pval, color = group)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Up (Dpysl3⁺)" = "red",
                                "Down (Dpysl3⁻)" = "blue",
                                "non-significant" = "gray")) +
  theme_bw(base_size = 14) +
  labs(title = "Volcano Plot: Cluster 5 Dpysl3⁺ vs Dpysl3⁻",
       x = "log2 Fold Change",
       y = "-log10(Adjusted P-value)",
       color = "Group") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  geom_text_repel(
    data = subset(deg, group != "non-significant"),
    aes(label = gene),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3
  )




# Gene list preparation for enrichment analyses

# Significance thresholds
padj_cutoff <- 0.05
logfc_cutoff <- 0.25

genes_up <- rownames(deg_cl5 %>%
                       filter(avg_log2FC > logfc_cutoff, p_val_adj < padj_cutoff))

length(genes_up)



# Gene Ontology (GO) enrichment analysis

# Convert gene symbols to Entrez IDs
gene_up_entrez <- bitr(
  genes_up,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)$ENTREZID


# GO enrichment for Cluster 5 cells
go_up_BP <- enrichGO(gene_up_entrez, org.Mm.eg.db, ont = "BP", readable = TRUE)
go_up_CC <- enrichGO(gene_up_entrez, org.Mm.eg.db, ont = "CC", readable = TRUE)
go_up_MF <- enrichGO(gene_up_entrez, org.Mm.eg.db, ont = "MF", readable = TRUE)

# KEGG pathway enrichment analysis
kegg_up <- enrichKEGG(gene_up_entrez, organism = "mmu")

dotplot(go_up_BP, showCategory = 15) +
  ggtitle("GO BP: Cluster 5")

dotplot(go_up_CC, showCategory = 15) +
  ggtitle("GO CC: Cluster 5")

dotplot(go_up_MF, showCategory = 15) +
  ggtitle("GO MF: Cluster 5")

dotplot(kegg_up, showCategory = 15) +
  ggtitle("KEGG: Cluster 5")

write.csv(as.data.frame(go_up_BP),
          "GO_BP_Cluster5.csv")

write.csv(as.data.frame(go_up_CC),
          "GO_CC_Cluster5.csv")

write.csv(as.data.frame(go_up_MF),
          "GO_MF_Cluster5.csv")

write.csv(as.data.frame(kegg_up),
          "KEGG_Cluster5.csv")



# Gene list preparation for enrichment analyses

# Significance thresholds

genes_up <- rownames(deg %>%
                       filter(avg_log2FC > logfc_cutoff, p_val_adj < padj_cutoff))

genes_down <- rownames(deg %>%
                         filter(avg_log2FC < -logfc_cutoff, p_val_adj < padj_cutoff))

length(genes_up)
length(genes_down)



# Gene Ontology (GO) enrichment analysis

# Convert gene symbols to Entrez IDs
gene_up_entrez <- bitr(
  genes_up,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)$ENTREZID

gene_down_entrez <- bitr(
  genes_down,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)$ENTREZID


# GO enrichment for Dpysl3-positive cells
go_up_BP <- enrichGO(gene_up_entrez, org.Mm.eg.db, ont = "BP", readable = TRUE)
go_up_CC <- enrichGO(gene_up_entrez, org.Mm.eg.db, ont = "CC", readable = TRUE)
go_up_MF <- enrichGO(gene_up_entrez, org.Mm.eg.db, ont = "MF", readable = TRUE)

# GO enrichment for Dpysl3-negative cells
go_down_BP <- enrichGO(gene_down_entrez, org.Mm.eg.db, ont = "BP", readable = TRUE)
go_down_CC <- enrichGO(gene_down_entrez, org.Mm.eg.db, ont = "CC", readable = TRUE)
go_down_MF <- enrichGO(gene_down_entrez, org.Mm.eg.db, ont = "MF", readable = TRUE)



# KEGG pathway enrichment analysis
kegg_up <- enrichKEGG(gene_up_entrez, organism = "mmu")
kegg_down <- enrichKEGG(gene_down_entrez, organism = "mmu")


dotplot(go_up_BP, showCategory = 15) +
  ggtitle("GO BP: Cluster5 Dpysl3⁺")

dotplot(go_down_BP, showCategory = 15) +
  ggtitle("GO BP: Cluster5 Dpysl3⁻")

dotplot(go_up_CC, showCategory = 15) +
  ggtitle("GO CC: Cluster5 Dpysl3⁺")

dotplot(go_down_CC, showCategory = 15) +
  ggtitle("GO CC: Cluster5 Dpysl3⁻")

dotplot(go_up_MF, showCategory = 15) +
  ggtitle("GO MF: Cluster5 Dpysl3⁺")

dotplot(go_down_MF, showCategory = 15) +
  ggtitle("GO MF: Cluster5 Dpysl3⁻")

dotplot(kegg_up, showCategory = 15) +
  ggtitle("KEGG: Cluster5 Dpysl3⁺")

dotplot(kegg_down, showCategory = 15) +
  ggtitle("KEGG: Cluster5 Dpysl3⁻")

write.csv(as.data.frame(go_up_BP),
          "GO_BP_Cluster5_Dpysl3_pos.csv")

write.csv(as.data.frame(go_down_BP),
          "GO_BP_Cluster5_Dpysl3_neg.csv")

write.csv(as.data.frame(go_up_CC),
          "GO_CC_Cluster5_Dpysl3_pos.csv")

write.csv(as.data.frame(go_down_CC),
          "GO_CC_Cluster5_Dpysl3_neg.csv")

write.csv(as.data.frame(go_up_MF),
          "GO_MF_Cluster5_Dpysl3_pos.csv")

write.csv(as.data.frame(go_down_MF),
          "GO_MF_Cluster5_Dpysl3_neg.csv")

write.csv(as.data.frame(kegg_up),
          "KEGG_Cluster5_Dpysl3_pos.csv")

write.csv(as.data.frame(kegg_down),
          "KEGG_Cluster5_Dpysl3_neg.csv")