#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Irene van Bloklad, Roy Oelen
# Name: epifat_read_batch1.R
# Function: read the first batch of spaceranger input into Seuart
############################################################################################################################


###############################
# loading libraries
###############################

library(Seurat) #"should be Seurat version 4!"
library(sctransform)
library(ggplot2)
library(spacexr)
library(Matrix)
library(cowplot)
library(dplyr)
library(patchwork)

###############################
# functions
###############################

get_spaceranger_objects <- function(space_ranger_loc, space_append="/outs/", spatial_image_loc="/spatial/") {
  # these are the sections
  sections <- list.dirs(space_ranger_loc, recursive = F, full.names = F)
  # get each of the sections
  sections_list <- list()
  for (section in sections) {
    # paste together the full directory
    space_full_dir <- paste(space_ranger_loc, section, space_append, sep = "")
    # and the image as well
    space_image_dir <- paste(space_full_dir, spatial_image_loc, sep = "")
    # and read the Visium image
    x10_image = Read10X_Image(space_image_dir)
    # the add the rest
    section_object <- Load10X_Spatial(space_full_dir, slice = section, image = x10_image)
    # add to the list
    sections_list[[section]] <- section_object
  }
  return(sections_list)
}

normalize_spaceranger_objects <- function(spaceranger_object_list) {
  # check each section
  for (section in names(spaceranger_object_list)) {
    # get the object from the list
    section_object <- spaceranger_object_list[[section]]
    # do whatever you want to do with the section object
    section_object <- SCTransform(section_object, assay = "Spatial")
    # put the object back in the list
    spaceranger_object_list[[section]] <- section_object
  }
  return(spaceranger_object_list)
}

make_violin_spaceranger_objects <- function(spaceranger_object_list) {
  # make a list to store the violins
  plot_list <- list()
  for (section in names(spaceranger_object_list)) {
    # get the object from the list
    section_object <- spaceranger_object_list[[section]]
    # set the identity so the plot is more clear
    section_object@meta.data[['section']] <- section
    # make the plot
    v_plot <- VlnPlot(section_object, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2, group.by = 'section') + NoLegend()
    # store the plot in the list
    plot_list[[section]] <- v_plot
  }
  p <- plot_grid(plotlist = plot_list)
  return(p)
}

filter_spaceranger_objects <- function(spaceranger_object_list) {
  # check each section
  for (section in names(spaceranger_object_list)) {
    # get the object from the list
    section_object <- spaceranger_object_list[[section]]
    # Filtering out  nFeature <100, nCount <200, and a gene expressed in <5 spots
    section_object <- section_object[, section_object$nFeature_Spatial > 100 & section_object$nCount_Spatial > 200 ]
    # put the object back in the list
    spaceranger_object_list[[section]] <- section_object
  }
  return(spaceranger_object_list)
}

prepare_spaceranger_objects_for_integration <- function(spaceranger_object_list) {
  # check each section
  for (section in names(spaceranger_object_list)) {
    # get the object from the list
    section_object <- spaceranger_object_list[[section]]
    # rename sections
    section_object@meta.data[["section"]] <- gsub("-", "_", section)
    section_object@meta.data[["barcode"]] <- colnames(section_object)
    section_object@meta.data[["section_barcode"]] <- paste(section_object@meta.data[["section"]], section_object@meta.data[["barcode"]], sep = "_")
    section_object <- RenameCells(section_object, new.names = section_object@meta.data[["section_barcode"]])
    
    # put the object back in the list
    spaceranger_object_list[[section]] <- section_object
  }
  return(spaceranger_object_list)
}

integrate_spaceranger_objects <- function(spaceranger_object_list) {
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

make_elbow_spaceranger_objects <- function(spaceranger_integrated) {
  # make a list to store the elbow
  plot_list <- list()
  # make the plot
  p <- ElbowPlot(object = spaceranger_integrated, ndims = 30)
  return(p)
}

make_umap_spaceranger_objects <- function(spaceranger_integrated) {
  # make a list to store the violins
  plot_list <- list()
  for (section in names(spaceranger_integrated)) {
    # get the object from the list
    section_object <- spaceranger_integrated[[section]]
    # set the identity so the plot is more clear
    section_object@meta.data[['section']] <- section
    # make the plot
    section_object <- RunUMAP(section_object, reduction = "pca", dims = 1:20)
    DimPlot(section_object, reduction = "umap", label = TRUE, label.size = 6,group.by = c("section"))
    # store the plot in the list
    plot_list[[section]] <- section_object
  }
  p <- plot_grid(plotlist = plot_list)
  return(p)
}

make_spatialumap_spaceranger_objects <- function(spaceranger_object_list) {
  # make a list to store the violins
  plot_list <- list()
  for (section in names(spaceranger_object_list)) {
    # get the object from the list
    section_object <- spaceranger_object_list[[section]]
    # set the identity so the plot is more clear
    section_object@meta.data[['section']] <- section
    # make the plot
    section_object <- SpatialDimPlot(section_object)
    # store the plot in the list
    plot_list[[section]] <- section_object
  }
  p <- plot_grid(plotlist = plot_list)
  return(p)
}


###############################
# order of data QC
###############################
# 1. Read data  - per section
# 2. Find UMIs and unique genes per section, make plot - per section
# 3. Filter out: A. Spots per section with too few genes, filter out genes if expressed in <5 spots?? - per section
# 4. Normalize data (SCTransform) - per section
# 5. intergrating data
# 6. PCA - merged
# 7. UMAP - merged
# 8. Deconvolution - per section using raw counts

###############################
# set the locations of the files
###############################
space_ranger_loc <- '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/epicardial_fat/processed/alignment/spaceranger_out/'

# location to store the resulting files
objects_loc <- "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/seurat_preprocess_samples/objects/"
spaceranger_object_loc <- paste(objects_loc, "spaceranger.20230811.rds", sep = "")

###############################
# 1. Read data - per section 
###############################
spaceranger_object_list <- get_spaceranger_objects(space_ranger_loc)

###############################
# 2. Find UMIs and unique genes per section, make plot - per section
###############################
make_violin_spaceranger_objects(spaceranger_object_list)

###############################
# 3. Filter out: nFeature_Spatial>100, nCount_Spatial>200 genes and ?expressed in <5 spots? - per section
###############################
spaceranger_object_list <- filter_spaceranger_objects(spaceranger_object_list)

###############################
# 4. Normalize data (SCTransform) - per section
###############################
spaceranger_object_list <- normalize_spaceranger_objects(spaceranger_object_list)
# save the object
saveRDS(spaceranger_object_list, spaceranger_object_loc)

###############################
# 5. merging and integrating data
###############################

# load newer Seurat version for integration
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
library(Seurat)

# prepare for integration
spaceranger_object_list <- prepare_spaceranger_objects_for_integration(spaceranger_object_list)
# integrating data
spaceranger_integrated <- integrate_spaceranger_objects(spaceranger_object_list)