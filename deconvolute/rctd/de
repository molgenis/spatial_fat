############################################################################################################################
# Authors: Irene van Blokland, Roy Oelen
# Name: epifat_de.R
# Function: differential expression analysis between cluster X vs all other
############################################################################################################################

###############################
# loading libraries
###############################
library(Seurat)

####################
# load files      
####################
# location of the file
objects_loc <- "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/seurat_preprocess_samples/objects/"
spaceranger_integrated_loc <- paste(objects_loc, "spaceranger.20230816.integrated.arteries.rds", sep = "")
# read the object
spaceranger_integrated <- readRDS(spaceranger_integrated_loc)
# output loc
de_output_loc <- "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/"

# make sure we have SCT normalized data
DefaultAssay(spaceranger_integrated) <- 'SCT'
# prepare DE 
spaceranger_integrated <- PrepSCTFindMarkers(spaceranger_integrated, assay = "SCT", verbose = TRUE)
# Do DE
markers <- FindAllMarkers(
  spaceranger_integrated,
  test.use = "wilcox",
  min.pct = 0.1,
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

# save DE lists as csv tables
write.csv(cluster_0_vs_all_others, file = "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/cluster_0_vs_all_others.csv", row.names = FALSE)
write.csv(cluster_1_vs_all_others, file = "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/cluster_1_vs_all_others.csv", row.names = FALSE)
write.csv(cluster_2_vs_all_others, file = "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/cluster_2_vs_all_others.csv", row.names = FALSE)
write.csv(cluster_3_vs_all_others, file = "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/cluster_3_vs_all_others.csv", row.names = FALSE)
write.csv(cluster_4_vs_all_others, file = "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/cluster_4_vs_all_others.csv", row.names = FALSE)
write.csv(cluster_5_vs_all_others, file = "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/cluster_5_vs_all_others.csv", row.names = FALSE)
write.csv(cluster_6_vs_all_others, file = "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/cluster_6_vs_all_others.csv", row.names = FALSE)
write.csv(cluster_7_vs_all_others, file = "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/cluster_7_vs_all_others.csv", row.names = FALSE)
write.csv(cluster_8_vs_all_others, file = "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/cluster_8_vs_all_others.csv", row.names = FALSE)
write.csv(cluster_9_vs_all_others, file = "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/cluster_9_vs_all_others.csv", row.names = FALSE)
write.csv(cluster_10_vs_all_others, file = "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/cluster_10_vs_all_others.csv", row.names = FALSE)
write.csv(cluster_11_vs_all_others, file = "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/cluster_11_vs_all_others.csv", row.names = FALSE)
write.csv(cluster_12_vs_all_others, file = "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/cluster_12_vs_all_others.csv", row.names = FALSE)
write.csv(cluster_13_vs_all_others, file = "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/cluster_13_vs_all_others.csv", row.names = FALSE)
write.csv(cluster_14_vs_all_others, file = "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/cluster_14_vs_all_others.csv", row.names = FALSE)

# downloading data to own computer
# readlink -f /groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/cluster_0_vs_all_others.csv
# on own computer in terminal
# scp tunnel+nibbler:/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/de/rctd/
  
