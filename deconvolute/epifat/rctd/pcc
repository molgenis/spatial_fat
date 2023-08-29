############################################################################################################################
# Authors: Irene van Blokland, Roy Oelen
# Name: epifat_analysis_template.R
# Function: perform PCC between predicted proportion of a specific cell type by RCTD and the expression profile of its marker gene
############################################################################################################################

###############################
# loading libraries
###############################

library(Seurat) #"should be Seurat version 5!"
library(Matrix)

###############################
# Calculating Pearsons´ correlation between predicted proportion of specific cell type and the expression profile of its marker gene
###############################

# read the objects
spaceranger_integrated <- readRDS("/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/seurat_preprocess_samples/objects/spaceranger.20230816.integrated.rds")
spaceranger_object_loc <- readRDS("/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/seurat_preprocess_samples/objects/spaceranger.20230811.rds")


############################################################################################################################
# PCC per cluster per section
############################################################################################################################

# count matrix selection of a specific section
sparse_matrix <- spaceranger_object_loc$`V10A20-016_A1`@assays$SCT@counts 
# make a matrix with spot x gene for a specific section 
dense_matrix <- as.matrix(sparse_matrix)
# Transpose the dense matrix to get spots x genes
spot_gene_matrix <- t(dense_matrix)
# Optionally, convert the transposed dense matrix back to a sparse matrix
spot_gene_sparse <- Matrix::Matrix(spot_gene_matrix, sparse = TRUE)
head(spot_gene_sparse)
## selecting my genes of interest for a specific section
genes_of_interest <- c("ADIPOQ", "PECAM1")
# Get the column indices corresponding to the genes of interest
gene_indices <- which(colnames(spot_gene_sparse) %in% genes_of_interest)
# Select columns with genes A and B
selected_matrix <- spot_gene_sparse[, gene_indices]

# fetch a few columns from the integrated RCTD deconvulated object
cell_type_prop <- spaceranger_integrated@meta.data[, c("Adipocyte", "Endothelial_cell", "seurat_clusters", "section", "barcode")]
# only continue with cells from cluster 0 and for a specific section
  # Filter rows based on conditions
filtered_rows <- cell_type_prop[cell_type_prop$seurat_clusters == 0 & cell_type_prop$section == "V10A20_016_A1", ]
  # Select specific columns
cell_prop_cluster_section <- filtered_rows[, c("Adipocyte", "Endothelial_cell", "barcode")]

# Extract spatial barcodes from the cell_type_prop
spatial_barcodes <- cell_prop_cluster_section$barcode
# Identify indices of spatial barcodes present in gene object
matching_indices <- which(rownames(selected_matrix) %in% spatial_barcodes)
# Subset the selected_matrix based on the matching indices
matched_sparse_matrix <- selected_matrix[matching_indices, ]
# Convert the matched_sparse_matrix to a dense matrix if needed
matched_dense_matrix <- as.matrix(matched_sparse_matrix)

# Merge to add columns from matched_dense_matrix to cell_prop_cluster_section
merged_all <- merge(cell_prop_cluster_section, matched_dense_matrix, by.x = "barcode", by.y = "row.names", all.x = TRUE)

# Calculate PCC Correlation according to benchmarking paper
gene_column <- merged_all$ADIPOQ
mean_gene_value <- mean(gene_column)
mean_cell_prop <- mean(merged_all$Adipocyte)
up <- sum((merged_all$ADIPOQ - mean_gene_value) * (merged_all$Adipocyte - mean_cell_prop))
down <- sqrt(sum((merged_all$ADIPOQ - mean_gene_value) ^ 2) * sum((merged_all$Adipocyte - mean_cell_prop) ^ 2))
rctd_corr <- up / down
cat("RCTD Corr:", rctd_corr, "\n")
# Adipocyte and ADIPOQ
corr <- cor.test(merged_all$Adipocyte, merged_all$ADIPOQ, method = "pearson")
print(corr)
# Endothelial_cell and PECAM1
corr <- cor.test(merged_all$Endothelial_cell, merged_all$PECAM1, method = "pearson")
print(corr)

## Pearsons correlation results per section per cluster and celltype
# Section A1, cluster 0
# Adipocyte and ADIPOQ,         corr = 0.3847437 , t = 4.5087, df = 117, p-value = 1.558e-05
# Endothelial cell and PECAM1,  corr = 0.02024948, t = 0.21908, df = 117, p-value = 0.827  

# Section B1, cluster 0
# Adipocyte and ADIPOQ,         corr = 0.3878969, t = 3.3668, df = 64, p-value = 0.001291
# Endothelial cell and PECAM1,  corr = 0.06915828, t = 0.55459, df = 64, p-value = 0.5811

# Section C1, cluster 0
# Adipocyte and ADIPOQ,         corr = 0.4380532, t = 9.8429, df = 408, p-value < 2.2e-16
# Endothelial cell and PECAM1,  corr = 0.2220115, t = 4.5992, df = 408, p-value = 5.668e-06

# Section D1, cluster 0
# Adipocyte and ADIPOQ,         corr = 0.3169707, t = 6.841, df = 419, p-value = 2.799e-11
# Endothelial cell and PECAM1,  corr = 0.003461646, t = 0.070859, df = 419, p-value = 0.9435




