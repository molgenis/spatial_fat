#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Irene van Bloklad, Roy Oelen
# Name: epifat_merge_batch1.R
# Function: merge batch1 sections
############################################################################################################################

###############################
# loading libraries
###############################

library(Seurat) #"should be Seurat version 5!"


###############################
# functions
###############################

integrate_spaceranger_objects <- function(spaceranger_object_list) {
  spaceranger_object_list <- lapply(X = spaceranger_object_list, FUN = SCTransform)
  # select integration features
  features <- SelectIntegrationFeatures(object.list = spaceranger_object_list, nfeatures = 3000)
  # create the anchors between the sections
  spaceranger_object_list <- PrepSCTIntegration(object.list = spaceranger_object_list, anchor.features = features)
  immune.anchors <- FindIntegrationAnchors(object.list = spaceranger_object_list, normalization.method = "SCT", anchor.features = features)
  # do integration
  se_merge_integrated <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
  # set the default assay
  DefaultAssay(se_merge_integrated) <- "integrated"
  return(spaceranger_object_list)
}


integrate_spaceranger_objects_rna_method <- function(spaceranger_object_list) {
  # normalize and identify variable features for each dataset independently
  spaceranger_object_list <- lapply(X = spaceranger_object_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = spaceranger_object_list)
  space_anchors <- FindIntegrationAnchors(object.list = spaceranger_object_list, anchor.features = features)
  space_combined <- IntegrateData(anchorset = space_anchors)
  # original unmodified data still resides in the 'Spatial' assay
  DefaultAssay(space_combined) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  space_combined <- ScaleData(space_combined, verbose = FALSE)
  space_combined <- RunPCA(space_combined, npcs = 30, verbose = FALSE)
  space_combined <- RunUMAP(space_combined, reduction = "pca", dims = 1:30)
  space_combined <- FindNeighbors(space_combined, reduction = "pca", dims = 1:30)
  space_combined <- FindClusters(space_combined, resolution = 1.2)
  return(space_combined)
}

###############################
# order of data QC
###############################
# 1. Read data  - per section
# 2. Find UMIs and unique genes per section, make plot - per section
# 3. Filter out: A. Spots per section with too few genes, filter out genes if expressed in <5 spots?? - per section
# 4. Normalize data (SCTransform) - per section
# 5. PCA/UMAP - merged
# 6. Deconvolution - per section using raw counts
# 7. merging and intergrating data

###############################
# set the locations of the files
###############################

# location to store the resulting files
objects_loc <- "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/seurat_preprocess_samples/objects/"
spaceranger_deconvolute_loc <- paste(objects_loc, "spaceranger.20230811.deconvoluted.rds", sep = "")
spaceranger_integrated_loc <- paste(objects_loc, "spaceranger.20230816.integrated.rds", sep = "")

###############################
# 6. merging and integrating data
###############################

# reda the object
spaceranger_object_list <- readRDS(spaceranger_deconvolute_loc)

# integrating all four sections
spaceranger_integrated <- integrate_spaceranger_objects_rna_method(spaceranger_object_list)

# add the deconvolution also a an assay
as_matrix <- t(spaceranger_integrated@meta.data[, 9:(ncol(spaceranger_integrated@meta.data) - 2)])
as_matrix <- matrix(as.numeric(as_matrix), ncol = ncol(as_matrix), dimnames = list(rownames(as_matrix), colnames(as_matrix)))
assay_object <- CreateAssay5Object(data = as_matrix)
spaceranger_integrated[['RCTD']] <- assay_object

# save integrated object
saveRDS(spaceranger_integrated, spaceranger_integrated_loc)
