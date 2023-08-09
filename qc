#############################
# Spatial transcriptomics data analysis 07-08-2023
# QC document of slice V10A20-016
# Using Seurat and ST QC pipeline https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html#Filter
#############################
# load packages and libraries
devtools::install_github("satijalab/seurat-data", dependencies = FALSE)
library(ggplot2)
library(dplyr)
library(Seurat)
#############################
# loading data and assembling everything into a table
# section 1
info_table <- data.frame(samples = c("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_A1/outs/filtered_feature_bc_matrix.h5"), 
                         imgs=c("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_A1/outs/spatial/tissue_lowres_image.png"), 
                         spotfiles=c("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_A1/outs/spatial/tissue_positions.csv"), 
                         json=c("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_A1/outs/spatial/scalefactors_json.json"))
# section 2
info_table_se2 <- data.frame(samples = c("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_B1/outs/filtered_feature_bc_matrix.h5"), 
                         imgs=c("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_B1/outs/spatial/tissue_lowres_image.png"), 
                         spotfiles=c("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_B1/outs/spatial/tissue_positions.csv"), 
                         json=c("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_B1/outs/spatial/scalefactors_json.json"))
# section 3
info_table_se3 <- data.frame(samples = c("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_C1/outs/filtered_feature_bc_matrix.h5"), 
                             imgs=c("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_C1/outs/spatial/tissue_lowres_image.png"), 
                             spotfiles=c("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_C1/outs/spatial/tissue_positions.csv"), 
                             json=c("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_C1/outs/spatial/scalefactors_json.json"))
# section 4
info_table_se4 <- data.frame(samples = c("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_D1/outs/filtered_feature_bc_matrix.h5"), 
                             imgs=c("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_D1/outs/spatial/tissue_lowres_image.png"), 
                             spotfiles=c("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_D1/outs/spatial/tissue_positions.csv"), 
                             json=c("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_D1/outs/spatial/scalefactors_json.json"))
#loading 10X Visium data into a Seurat object
se_seurat <- Load10X_Spatial("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_A1/outs/", image = Read10X_Image("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_A1/outs/spatial/"), slice = "A1")
se2_seurat <- Load10X_Spatial("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_B1/outs/", slice = "B1")
se3_seurat <- Load10X_Spatial("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_C1/outs/", slice = "C1")
se4_seurat <- Load10X_Spatial("/Volumes/RoyExtSSD1T/work/genetica/epicardial_fat/processed/alignment/spaceranger_out/V10A20-016_D1/outs/", slice = "D1")
#############################
# Step 1: plotting number of genes, Hb and mito genes: These are not present in the ST probe set and so not calculated.
# generating violin plots to see number of molecules and genes
vl1 <- VlnPlot(se_seurat, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
vl2 <- VlnPlot(se2_seurat, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
vl3 <- VlnPlot(se3_seurat, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
vl4 <- VlnPlot(se4_seurat, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
vl1+vl2+vl3+vl4
#############################
# Step 2: Filtering out spots with too few genes
se_seurat <- se_seurat[, se_seurat$nFeature_Spatial > 100] # for MT/Hb filter add: & se$percent_mito < 25 & brain$percent_hb <20]
se2_seurat <- se2_seurat[, se2_seurat$nFeature_Spatial > 100] 
se3_seurat <- se3_seurat[, se3_seurat$nFeature_Spatial > 100] 
se4_seurat <- se4_seurat[, se4_seurat$nFeature_Spatial > 100] 
#############################
# Step 3: merging all 4 sections into the dataset
se_seurat@meta.data[["section"]] <- "A1"
se2_seurat@meta.data[["section"]] <- "B1"
se3_seurat@meta.data[["section"]] <- "C1"
se4_seurat@meta.data[["section"]] <- "D1"
se_seurat@meta.data[["barcode"]] <- colnames(se_seurat)
se2_seurat@meta.data[["barcode"]] <- colnames(se2_seurat)
se3_seurat@meta.data[["barcode"]] <- colnames(se3_seurat)
se4_seurat@meta.data[["barcode"]] <- colnames(se4_seurat)
se_seurat@meta.data[["section_barcode"]] <- paste(se_seurat@meta.data[["section"]], se_seurat@meta.data[["barcode"]], sep = "_")
se2_seurat@meta.data[["section_barcode"]] <- paste(se2_seurat@meta.data[["section"]], se2_seurat@meta.data[["barcode"]], sep = "_")
se3_seurat@meta.data[["section_barcode"]] <- paste(se3_seurat@meta.data[["section"]], se3_seurat@meta.data[["barcode"]], sep = "_")
se4_seurat@meta.data[["section_barcode"]] <- paste(se4_seurat@meta.data[["section"]], se4_seurat@meta.data[["barcode"]], sep = "_")
se_seurat <- RenameCells(se_seurat, new.names = se_seurat@meta.data[["section_barcode"]])
se2_seurat <- RenameCells(se2_seurat, new.names = se2_seurat@meta.data[["section_barcode"]])
se3_seurat <- RenameCells(se3_seurat, new.names = se3_seurat@meta.data[["section_barcode"]])
se4_seurat <- RenameCells(se4_seurat, new.names = se4_seurat@meta.data[["section_barcode"]])
se_merge <- merge(se_seurat, c(se2_seurat, se3_seurat, se4_seurat))
#############################
# Step 4: integrating the data
st.list = list(A1 = se_seurat, B1= se2_seurat, C1= se3_seurat, D1=se4_seurat)
features <- SelectIntegrationFeatures(object.list = st.list, verbose = FALSE)
immune.anchors <- FindIntegrationAnchors(object.list = st.list, anchor.features = features)
# this command creates an ‘integrated’ data assay
se_merge_integrated <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(se_merge_integrated) <- "integrated"
#############################
# Step 5: checking for PCs
ElbowPlot(object = se_merge_integrated, ndims = 30)
#############################
# Step 6: Dimensionality reduction and clustering
se_merge_integrated <- ScaleData(se_merge_integrated, verbose = FALSE)
se_merge_integrated <- RunPCA(se_merge_integrated, npcs = 30, verbose = FALSE)
se_merge_integrated <- FindNeighbors(se_merge_integrated, reduction = "pca", dims = 1:20)
se_merge_integrated <- FindClusters(se_merge_integrated, resolution = 1.2)
se_merge_integrated <- RunUMAP(se_merge_integrated, reduction = "pca", dims = 1:20)
p1 <- DimPlot(se_merge_integrated, reduction = "umap")
p2 <- DimPlot(se_merge_integrated, reduction = "umap", label = TRUE, repel = TRUE)
# plot the DimPlots
p1 + p2
# plot the clusters on the H&E image
SpatialDimPlot(se_merge_integrated)
#############################
# Step 7: plotting specific features to see whether clustering makes sense
# normalize the actual genes as well
se_merge_integrated <- NormalizeData(se_merge_integrated, assay = "Spatial")
# set that as the active assay
DefaultAssay(se_merge_integrated) <- "Spatial"
SpatialFeaturePlot(se_merge_integrated, features=c("IL6"))
SpatialFeaturePlot(se_merge_integrated, features=c("CD27"))
plot_grid(SpatialFeaturePlot(se_merge_integrated, features=c("ADIPOQ")), SpatialDimPlot(se_merge_integrated, group.by = "seurat_clusters"), plot_grid(FeaturePlot(se_merge_integrated, features = c("ADIPOQ")), DimPlot(se_merge_integrated, group.by = "seurat_clusters"), ncol=2), nrow = 3)
#############################
# Step 8: Deconvolution: Robust decomposition of cell type mixtures in spatial transcriptomics 
