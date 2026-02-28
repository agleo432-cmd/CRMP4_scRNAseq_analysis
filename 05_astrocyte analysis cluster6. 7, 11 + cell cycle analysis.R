# This script identifies CRMP4 (Dpysl3)–positive and –negative cells within 
# astrocyte clusters 6 and 11, performs differential gene expression 
# analysis between these two states, and visualizes the results using a 
# volcano plot. Gene Ontology (BP, CC, MF) and KEGG pathway enrichment analyses
# are then conducted to characterize functional differences associated with 
# Dpysl3 expression.Cell cycle phase scoring was additionally done and
# was overlaped on the UMAP.

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

# for loading
astro_merged <- readRDS("integrated_data_astro_injured_GSE189070.rds")

#########   Clustering    ###########################################

DefaultAssay(astro_merged) <- "integrated"

astro_merged <- ScaleData(astro_merged)
astro_merged <- RunPCA(astro_merged)
astro_merged <- RunUMAP(astro_merged, dims = 1:20)

astro_merged <- FindNeighbors(astro_merged,dims=1:30)
astro_merged <- FindClusters(astro_merged, resolution = c(0.4))

DimPlot(astro_merged, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP (Resolution 0.4)") +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed()



##### Expression & function analysis of CRMP4 expressing astrocytes (Cluster 6)  #########

############################################################
# Analysis workflow:
# Differential expression and functional enrichment analysis
# based on CRMP4 (Dpysl3) expression status
############################################################

# Define CRMP4 (Dpysl3) expression status in subset clusters

# extract cluster 6 astrocytes
cluster6 <- subset(astro_merged, idents = "6")

gene_name <- "Dpysl3"
expr_mat <- GetAssayData(cluster6, slot = "data")

Dpysl3_pos <- expr_mat[gene_name, ] > 0

cluster6$Dpysl3_status <- ifelse(Dpysl3_pos, "Dpysl3_pos", "Dpysl3_neg")

table(cluster6$Dpysl3_status)

# Differential gene expression analysis (Cluster 6)

deg_cl6 <- FindMarkers(
  astro_merged,
  ident.1 = 6,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

head(deg_cl6)

# save as csv
write.csv(deg_cl6, file = "DEG cluster6 hemi SCI astrocyte.csv")

deg_cl6$log10_pval_plot <- pmin(
  -log10(deg_cl6$p_val_adj),
  300
)

# Annotate DEG directionality
deg_cl6$group <- "non-significant"
deg_cl6$group[deg_cl6$avg_log2FC > 0.25 & deg_cl6$p_val_adj < 0.05] <- "Up (Cluster6)"
deg_cl6$group[deg_cl6$avg_log2FC < -0.25 & deg_cl6$p_val_adj < 0.05] <- "Down"

deg_cl6$gene <- rownames(deg_cl6)

ggplot(deg_cl6, aes(x = avg_log2FC, y = log10_pval_plot, color = group)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c(
    "Up (Cluster6)" = "red",
    "Down" = "blue",
    "non-significant" = "gray"
  )) +
  theme_bw(base_size = 14) +
  labs(
    title = "Volcano Plot: Cluster 6",
    x = "log2 Fold Change",
    y = "-log10(Adjusted P-value)",
    color = "Group"
  ) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = subset(deg_cl6, group != "non-significant"),
    aes(label = gene),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3
  )

#  Differential gene expression analysis (Dpysl3+ vs Dpysl3− in cluster 6)

Idents(cluster6) <- "Dpysl3_status"

# Detecting DEGs
Dpysl3_DEG <- FindMarkers(cluster6,
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

# save as csv
write.csv(deg, file = "DEG_cluster6 hemi SCI astrocyte CRMP4+.csv")

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
  labs(title = "Volcano Plot: Cluster 6 Dpysl3⁺ vs Dpysl3⁻",
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

genes_up <- rownames(deg_cl6 %>%
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


# GO enrichment for Cluster 6 cells
go_up_BP <- enrichGO(gene_up_entrez, org.Mm.eg.db, ont = "BP", readable = TRUE)
go_up_CC <- enrichGO(gene_up_entrez, org.Mm.eg.db, ont = "CC", readable = TRUE)
go_up_MF <- enrichGO(gene_up_entrez, org.Mm.eg.db, ont = "MF", readable = TRUE)

# KEGG pathway enrichment analysis
kegg_up <- enrichKEGG(gene_up_entrez, organism = "mmu")

dotplot(go_up_BP, showCategory = 15) +
  ggtitle("GO BP: Cluster 6")

dotplot(go_up_CC, showCategory = 15) +
  ggtitle("GO CC: Cluster 6")

dotplot(go_up_MF, showCategory = 15) +
  ggtitle("GO MF: Cluster 6")

dotplot(kegg_up, showCategory = 15) +
  ggtitle("KEGG: Cluster 6")

write.csv(as.data.frame(go_up_BP),
          "GO_BP_Cluster6.csv")

write.csv(as.data.frame(go_up_CC),
          "GO_CC_Cluster6.csv")

write.csv(as.data.frame(go_up_MF),
          "GO_MF_Cluster6.csv")

write.csv(as.data.frame(kegg_up),
          "KEGG_Cluster6.csv")



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
  ggtitle("GO BP: Cluster6 Dpysl3⁺")

dotplot(go_down_BP, showCategory = 15) +
  ggtitle("GO BP: Cluster6 Dpysl3⁻")

dotplot(go_up_CC, showCategory = 15) +
  ggtitle("GO CC: Cluster6 Dpysl3⁺")

dotplot(go_down_CC, showCategory = 15) +
  ggtitle("GO CC: Cluster6 Dpysl3⁻")

dotplot(go_up_MF, showCategory = 15) +
  ggtitle("GO MF: Cluster6 Dpysl3⁺")

dotplot(go_down_MF, showCategory = 15) +
  ggtitle("GO MF: Cluster6 Dpysl3⁻")

dotplot(kegg_up, showCategory = 15) +
  ggtitle("KEGG: Cluster6 Dpysl3⁺")

dotplot(kegg_down, showCategory = 15) +
  ggtitle("KEGG: Cluster6 Dpysl3⁻")

write.csv(as.data.frame(go_up_BP),
          "GO_BP_Cluster6_Dpysl3_pos.csv")

write.csv(as.data.frame(go_down_BP),
          "GO_BP_Cluster6_Dpysl3_neg.csv")

write.csv(as.data.frame(go_up_CC),
          "GO_CC_Cluster6_Dpysl3_pos.csv")

write.csv(as.data.frame(go_down_CC),
          "GO_CC_Cluster6_Dpysl3_neg.csv")

write.csv(as.data.frame(go_up_MF),
          "GO_MF_Cluster6_Dpysl3_pos.csv")

write.csv(as.data.frame(go_down_MF),
          "GO_MF_Cluster6_Dpysl3_neg.csv")

write.csv(as.data.frame(kegg_up),
          "KEGG_Cluster6_Dpysl3_pos.csv")

write.csv(as.data.frame(kegg_down),
          "KEGG_Cluster6_Dpysl3_neg.csv")


###############################################################################

# 1. Define CRMP4 (Dpysl3) expression status in subset clusters

# extract cluster 7 astrocytes
cluster7 <- subset(astro_merged, idents = "7")

gene_name <- "Dpysl3"
expr_mat <- GetAssayData(cluster7, slot = "data")

Dpysl3_pos <- expr_mat[gene_name, ] > 0

cluster7$Dpysl3_status <- ifelse(Dpysl3_pos, "Dpysl3_pos", "Dpysl3_neg")

table(cluster7$Dpysl3_status)

# 2. Differential gene expression analysis (Cluster 7)

deg_cl7 <- FindMarkers(
  astro_merged,
  ident.1 = 7,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

head(deg_cl7)

# save as csv
write.csv(deg_cl7, file = "DEG cluster7 hemi SCI astrocyte.csv")

deg_cl7$log10_pval_plot <- pmin(
  -log10(deg_cl7$p_val_adj),
  300
)

# Annotate DEG directionality
deg_cl7$group <- "non-significant"
deg_cl7$group[deg_cl7$avg_log2FC > 0.25 & deg_cl7$p_val_adj < 0.05] <- "Up (Cluster7)"
deg_cl7$group[deg_cl7$avg_log2FC < -0.25 & deg_cl7$p_val_adj < 0.05] <- "Down"

deg_cl7$gene <- rownames(deg_cl7)

ggplot(deg_cl7, aes(x = avg_log2FC, y = log10_pval_plot, color = group)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c(
    "Up (Cluster7)" = "red",
    "Down" = "blue",
    "non-significant" = "gray"
  )) +
  theme_bw(base_size = 14) +
  labs(
    title = "Volcano Plot: Cluster 7",
    x = "log2 Fold Change",
    y = "-log10(Adjusted P-value)",
    color = "Group"
  ) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = subset(deg_cl7, group != "non-significant"),
    aes(label = gene),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3
  )

# Differential gene expression analysis (Dpysl3+ vs Dpysl3−)

Idents(cluster7) <- "Dpysl3_status"

# Detecting DEGs
Dpysl3_DEG <- FindMarkers(cluster7,
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

# save as csv
write.csv(deg, file = "DEG_cluster7 hemi SCI astrocyte CRMP4+.csv")

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
  labs(title = "Volcano Plot: Cluster 7 Dpysl3⁺ vs Dpysl3⁻",
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

genes_up <- rownames(deg_cl7 %>%
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
  ggtitle("GO BP: Cluster 7")

dotplot(go_up_CC, showCategory = 15) +
  ggtitle("GO CC: Cluster 7")

dotplot(go_up_MF, showCategory = 15) +
  ggtitle("GO MF: Cluster 7")

dotplot(kegg_up, showCategory = 15) +
  ggtitle("KEGG: Cluster 7")

write.csv(as.data.frame(go_up_BP),
          "GO_BP_Cluster7.csv")

write.csv(as.data.frame(go_up_CC),
          "GO_CC_Cluster7.csv")

write.csv(as.data.frame(go_up_MF),
          "GO_MF_Cluster7.csv")

write.csv(as.data.frame(kegg_up),
          "KEGG_Cluster7.csv")



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
  ggtitle("GO BP: Cluster7 Dpysl3⁺")

dotplot(go_down_BP, showCategory = 15) +
  ggtitle("GO BP: Cluster7 Dpysl3⁻")

dotplot(go_up_CC, showCategory = 15) +
  ggtitle("GO CC: Cluster7 Dpysl3⁺")

dotplot(go_down_CC, showCategory = 15) +
  ggtitle("GO CC: Cluster7 Dpysl3⁻")

dotplot(go_up_MF, showCategory = 15) +
  ggtitle("GO MF: Cluster7 Dpysl3⁺")

dotplot(go_down_MF, showCategory = 15) +
  ggtitle("GO MF: Cluster7 Dpysl3⁻")

dotplot(kegg_up, showCategory = 15) +
  ggtitle("KEGG: Cluster7 Dpysl3⁺")

dotplot(kegg_down, showCategory = 15) +
  ggtitle("KEGG: Cluster7 Dpysl3⁻")

write.csv(as.data.frame(go_up_BP),
          "GO_BP_Cluster7_Dpysl3_pos.csv")

write.csv(as.data.frame(go_down_BP),
          "GO_BP_Cluster7_Dpysl3_neg.csv")

write.csv(as.data.frame(go_up_CC),
          "GO_CC_Cluster7_Dpysl3_pos.csv")

write.csv(as.data.frame(go_down_CC),
          "GO_CC_Cluster7_Dpysl3_neg.csv")

write.csv(as.data.frame(go_up_MF),
          "GO_MF_Cluster7_Dpysl3_pos.csv")

write.csv(as.data.frame(go_down_MF),
          "GO_MF_Cluster7_Dpysl3_neg.csv")

write.csv(as.data.frame(kegg_up),
          "KEGG_Cluster7_Dpysl3_pos.csv")

write.csv(as.data.frame(kegg_down),
          "KEGG_Cluster7_Dpysl3_neg.csv")


###############################################################################

# Define CRMP4 (Dpysl3) expression status in subset clusters

# extract cluster 11 astrocytes
cluster11 <- subset(astro_merged, idents = "11")

gene_name <- "Dpysl3"
expr_mat <- GetAssayData(cluster11, slot = "data")

Dpysl3_pos <- expr_mat[gene_name, ] > 0

cluster11$Dpysl3_status <- ifelse(Dpysl3_pos, "Dpysl3_pos", "Dpysl3_neg")

table(cluster11$Dpysl3_status)

# Differential gene expression analysis (Cluster 11)

deg_cl11 <- FindMarkers(
  astro_merged,
  ident.1 = 11,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

head(deg_cl11)

# save as csv
write.csv(deg_cl11, file = "DEG cluster11 hemi SCI astrocyte.csv")

deg_cl11$log10_pval_plot <- pmin(
  -log10(deg_cl11$p_val_adj),
  300
)

# Annotate DEG directionality
deg_cl11$group <- "non-significant"
deg_cl11$group[deg_cl11$avg_log2FC > 0.25 & deg_cl11$p_val_adj < 0.05] <- "Up (Cluster11)"
deg_cl11$group[deg_cl11$avg_log2FC < -0.25 & deg_cl11$p_val_adj < 0.05] <- "Down"

deg_cl11$gene <- rownames(deg_cl11)

ggplot(deg_cl11, aes(x = avg_log2FC, y = log10_pval_plot, color = group)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c(
    "Up (Cluster11)" = "red",
    "Down" = "blue",
    "non-significant" = "gray"
  )) +
  theme_bw(base_size = 14) +
  labs(
    title = "Volcano Plot: Cluster 11",
    x = "log2 Fold Change",
    y = "-log10(Adjusted P-value)",
    color = "Group"
  ) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = subset(deg_cl11, group != "non-significant"),
    aes(label = gene),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3
  )

# Differential gene expression analysis (Dpysl3+ vs Dpysl3−)

Idents(cluster11) <- "Dpysl3_status"

# Detecting DEGs
Dpysl3_DEG <- FindMarkers(cluster11,
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

# save as csv
write.csv(deg, file = "DEG_cluster11 hemi SCI astrocyte CRMP4+.csv")

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
  labs(title = "Volcano Plot: Cluster 11 Dpysl3⁺ vs Dpysl3⁻",
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

genes_up <- rownames(deg_cl11 %>%
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
  ggtitle("GO BP: Cluster 11")

dotplot(go_up_CC, showCategory = 15) +
  ggtitle("GO CC: Cluster 11")

dotplot(go_up_MF, showCategory = 15) +
  ggtitle("GO MF: Cluster 11")

dotplot(kegg_up, showCategory = 15) +
  ggtitle("KEGG: Cluster 11")

write.csv(as.data.frame(go_up_BP),
          "GO_BP_Cluster11.csv")

write.csv(as.data.frame(go_up_CC),
          "GO_CC_Cluster11.csv")

write.csv(as.data.frame(go_up_MF),
          "GO_MF_Cluster11.csv")

write.csv(as.data.frame(kegg_up),
          "KEGG_Cluster11.csv")



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
  ggtitle("GO BP: Cluster11 Dpysl3⁺")

dotplot(go_down_BP, showCategory = 15) +
  ggtitle("GO BP: Cluster11 Dpysl3⁻")

dotplot(go_up_CC, showCategory = 15) +
  ggtitle("GO CC: Cluster11 Dpysl3⁺")

dotplot(go_down_CC, showCategory = 15) +
  ggtitle("GO CC: Cluster11 Dpysl3⁻")

dotplot(go_up_MF, showCategory = 15) +
  ggtitle("GO MF: Cluster11 Dpysl3⁺")

dotplot(go_down_MF, showCategory = 15) +
  ggtitle("GO MF: Cluster11 Dpysl3⁻")

dotplot(kegg_up, showCategory = 15) +
  ggtitle("KEGG: Cluster11 Dpysl3⁺")

dotplot(kegg_down, showCategory = 15) +
  ggtitle("KEGG: Cluster11 Dpysl3⁻")

write.csv(as.data.frame(go_up_BP),
          "GO_BP_Cluster11_Dpysl3_pos.csv")

write.csv(as.data.frame(go_down_BP),
          "GO_BP_Cluster11_Dpysl3_neg.csv")

write.csv(as.data.frame(go_up_CC),
          "GO_CC_Cluster11_Dpysl3_pos.csv")

write.csv(as.data.frame(go_down_CC),
          "GO_CC_Cluster11_Dpysl3_neg.csv")

write.csv(as.data.frame(go_up_MF),
          "GO_MF_Cluster11_Dpysl3_pos.csv")

write.csv(as.data.frame(go_down_MF),
          "GO_MF_Cluster11_Dpysl3_neg.csv")

write.csv(as.data.frame(kegg_up),
          "KEGG_Cluster11_Dpysl3_pos.csv")

write.csv(as.data.frame(kegg_down),
          "KEGG_Cluster11_Dpysl3_neg.csv")


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

p1
