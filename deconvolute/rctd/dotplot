####################
# make clustered dotplot
####################
library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork) 
library(Seurat)

# prepare DE 
spaceranger_object <- PrepSCTFindMarkers(spaceranger_object, assay = "SCT", verbose = TRUE)

# Do DE
markers <- FindAllMarkers(
  spaceranger_object,
  test.use = "wilcox",
  min.pct = 0.1,
  assay = "SCT"
)

# Fetch the markers lists for every cluster versus all other clusters
cluster_0_vs_all_others <- markers[markers$cluster == 0, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
cluster_1_vs_all_others <- markers[markers$cluster == 1, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
cluster_2_vs_all_others <- markers[markers$cluster == 2, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
cluster_3_vs_all_others <- markers[markers$cluster == 3, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
cluster_4_vs_all_others <- markers[markers$cluster == 4, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
cluster_5_vs_all_others <- markers[markers$cluster == 5, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
cluster_6_vs_all_others <- markers[markers$cluster == 6, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
cluster_7_vs_all_others <- markers[markers$cluster == 7, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
cluster_8_vs_all_others <- markers[markers$cluster == 8, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
cluster_9_vs_all_others <- markers[markers$cluster == 9, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
cluster_10_vs_all_others <- markers[markers$cluster == 10, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
cluster_11_vs_all_others <- markers[markers$cluster == 11, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
cluster_12_vs_all_others <- markers[markers$cluster == 12, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
cluster_13_vs_all_others <- markers[markers$cluster == 13, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
cluster_14_vs_all_others <- markers[markers$cluster == 14, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]


# Sort the DEGs by adjusted p-value (assuming it's in the 'p_value' column)
cluster_0_vs_all_others_sorted <- cluster_0_vs_all_others[order(cluster_0_vs_all_others$p_val_adj), ]
cluster_1_vs_all_others_sorted <- cluster_1_vs_all_others[order(cluster_1_vs_all_others$p_val_adj), ]
cluster_2_vs_all_others_sorted <- cluster_2_vs_all_others[order(cluster_2_vs_all_others$p_val_adj), ]
cluster_3_vs_all_others_sorted <- cluster_3_vs_all_others[order(cluster_3_vs_all_others$p_val_adj), ]
cluster_4_vs_all_others_sorted <- cluster_4_vs_all_others[order(cluster_4_vs_all_others$p_val_adj), ]
cluster_5_vs_all_others_sorted <- cluster_5_vs_all_others[order(cluster_5_vs_all_others$p_val_adj), ]
cluster_6_vs_all_others_sorted <- cluster_6_vs_all_others[order(cluster_6_vs_all_others$p_val_adj), ]
cluster_7_vs_all_others_sorted <- cluster_7_vs_all_others[order(cluster_7_vs_all_others$p_val_adj), ]
cluster_8_vs_all_others_sorted <- cluster_8_vs_all_others[order(cluster_8_vs_all_others$p_val_adj), ]
cluster_9_vs_all_others_sorted <- cluster_9_vs_all_others[order(cluster_9_vs_all_others$p_val_adj), ]
cluster_10_vs_all_others_sorted <- cluster_10_vs_all_others[order(cluster_10_vs_all_others$p_val_adj), ]
cluster_11_vs_all_others_sorted <- cluster_11_vs_all_others[order(cluster_11_vs_all_others$p_val_adj), ]
cluster_12_vs_all_others_sorted <- cluster_12_vs_all_others[order(cluster_12_vs_all_others$p_val_adj), ]
cluster_13_vs_all_others_sorted <- cluster_13_vs_all_others[order(cluster_13_vs_all_others$p_val_adj), ]
cluster_14_vs_all_others_sorted <- cluster_14_vs_all_others[order(cluster_14_vs_all_others$p_val_adj), ]

# Fetch the top 5 DEGs and add a cluster column
top_5_genes_cluster0 <- head(cluster_0_vs_all_others_sorted, n = 5)
top_5_genes_cluster0$cluster <- 0
top_5_genes_cluster1 <- head(cluster_1_vs_all_others_sorted, n = 5)
top_5_genes_cluster1$cluster <- 1
top_5_genes_cluster2 <- head(cluster_2_vs_all_others_sorted, n = 5)
top_5_genes_cluster2$cluster <- 2
top_5_genes_cluster3 <- head(cluster_3_vs_all_others_sorted, n = 5)
top_5_genes_cluster3$cluster <- 3
top_5_genes_cluster4 <- head(cluster_4_vs_all_others_sorted, n = 5)
top_5_genes_cluster4$cluster <- 4
top_5_genes_cluster5 <- head(cluster_5_vs_all_others_sorted, n = 5)
top_5_genes_cluster5$cluster <- 5
top_5_genes_cluster6 <- head(cluster_6_vs_all_others_sorted, n = 5)
top_5_genes_cluster6$cluster <- 6
top_5_genes_cluster7 <- head(cluster_7_vs_all_others_sorted, n = 5)
top_5_genes_cluster7$cluster <- 7
top_5_genes_cluster8 <- head(cluster_8_vs_all_others_sorted, n = 5)
top_5_genes_cluster8$cluster <- 8
top_5_genes_cluster9 <- head(cluster_9_vs_all_others_sorted, n = 5)
top_5_genes_cluster9$cluster <- 9
top_5_genes_cluster10 <- head(cluster_10_vs_all_others_sorted, n = 5)
top_5_genes_cluster10$cluster <- 10
top_5_genes_cluster11 <- head(cluster_11_vs_all_others_sorted, n = 5)
top_5_genes_cluster11$cluster <- 11
top_5_genes_cluster12 <- head(cluster_12_vs_all_others_sorted, n = 5)
top_5_genes_cluster12$cluster <- 12
top_5_genes_cluster13 <- head(cluster_13_vs_all_others_sorted, n = 5)
top_5_genes_cluster13$cluster <- 13
top_5_genes_cluster14 <- head(cluster_14_vs_all_others_sorted, n = 5)
top_5_genes_cluster14$cluster <- 14

# merging all top genes
merged_top_genes <- rbind(top_5_genes_cluster0, top_5_genes_cluster1, top_5_genes_cluster2,
                          top_5_genes_cluster3, top_5_genes_cluster4, top_5_genes_cluster5, 
                          top_5_genes_cluster6, top_5_genes_cluster7, top_5_genes_cluster8, 
                          top_5_genes_cluster9, top_5_genes_cluster10, top_5_genes_cluster11,
                          top_5_genes_cluster12, top_5_genes_cluster13, top_5_genes_cluster14)

# extract gene names
genes <- merged_top_genes$gene

# prepare data
cluster_0_vs_all_others_sorted$cluster <- 0
cluster_1_vs_all_others_sorted$cluster <- 1
cluster_2_vs_all_others_sorted$cluster <- 2
cluster_3_vs_all_others_sorted$cluster <- 3
cluster_4_vs_all_others_sorted$cluster <- 4
cluster_5_vs_all_others_sorted$cluster <- 5
cluster_6_vs_all_others_sorted$cluster <- 6
cluster_7_vs_all_others_sorted$cluster <- 7
cluster_8_vs_all_others_sorted$cluster <- 8
cluster_9_vs_all_others_sorted$cluster <- 9
cluster_10_vs_all_others_sorted$cluster <- 10
cluster_11_vs_all_others_sorted$cluster <- 11
cluster_12_vs_all_others_sorted$cluster <- 12
cluster_13_vs_all_others_sorted$cluster <- 13
cluster_14_vs_all_others_sorted$cluster <- 14

# merging
merged_genes <- rbind(cluster_0_vs_all_others_sorted, cluster_1_vs_all_others_sorted, cluster_2_vs_all_others_sorted,
                      cluster_3_vs_all_others_sorted, cluster_4_vs_all_others_sorted, cluster_5_vs_all_others_sorted,
                      cluster_6_vs_all_others_sorted, cluster_7_vs_all_others_sorted, cluster_8_vs_all_others_sorted,
                      cluster_9_vs_all_others_sorted, cluster_10_vs_all_others_sorted, cluster_11_vs_all_others_sorted,
                      cluster_12_vs_all_others_sorted, cluster_13_vs_all_others_sorted, cluster_14_vs_all_others_sorted)

# vector of top 5 DE genes per cluster
#top_genes_list <- data.frame(merged_top_genes$gene)
# extract rows for all genes that were in the top five at some point
merged_genes_tops <- merged_genes[merged_genes$gene %in% genes, ]
# remove duplicates
genes_nondup <- genes[!duplicated(genes)]
# make the genes a factor of our given order
merged_genes_tops$gene <- factor(merged_genes_tops$gene, levels = genes_nondup)
# make the zero p values very small
merged_genes_tops$pval_safe <- merged_genes_tops$p_val_adj
merged_genes_tops[merged_genes_tops$pval_safe == 0, "pval_safe"] <- .Machine$double.xmin
# log10 transform
merged_genes_tops$pval_safe_log <- -log10(merged_genes_tops$pval_safe)
# make the numbers characters, so they are all shown
merged_genes_tops$cluster <- factor(as.character(merged_genes_tops$cluster), levels = unique(as.character(merged_genes_tops$cluster)))

# make the plot
ggplot(merged_genes_tops, aes(x = cluster, y = gene, size = pval_safe_log, color = avg_log2FC)) +
  geom_point() +
  scale_colour_gradient2(
    high = "darkred",
    mid = "grey",
    low = "darkblue",
    midpoint = 0,
    guide = guide_colorbar(nbin = 4),
    limits = c(-5.2,5.2),
) + scale_size_continuous(range = c(1, 4)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.border = element_rect(color="black", fill=NA, size=1.1), 
          panel.grid.major = element_line(color = "gray", size = 0.25, linetype = 1), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          strip.background = element_rect(colour="white", fill="white")) +
  labs(size = "-log10(p)")
