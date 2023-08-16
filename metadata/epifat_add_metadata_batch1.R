#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: epifat_add_metadata_batch1.R
# Function: add metadata for each section
############################################################################################################################

###############################
# loading libraries
###############################

library(Seurat) #"should be Seurat version 5!"


###############################
# functions
###############################

#' add metadata to a Seurat object
#' 
#' @param seurat_object the Seurat object to add the metadata to
#' @param section_column the column that denotes the section of the spot
#' @param barcode_column the column that denotes the 10x barcode of the spot
#' @param interested_metadata_columns which metadata columns to add
#' @returns the Seurat object with the metadata added
#' 
add_metadata <- function(seurat_object, metadata_table, section_column='section', barcode_column='barcode', interested_metadata_columns=c('artery', 'status')) {
  # add rownames to the metadata table, that reflect the combination of section and barcode
  rownames(metadata_table) <- paste(metadata_table[[section_column]], metadata_table[[barcode_column]], sep = '_')
  # next, subset to only the columns we care about
  metadata_table <- metadata_table[, interested_metadata_columns]
  # now add it to the Seurat object
  seurat_object <- AddMetaData(seurat_object, metadata_table)
  # and return the result
  return(seurat_object)
}


###############################
# settings
###############################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

###############################
# set the locations of the files
###############################

# location to store the resulting files
objects_loc <- "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/seurat_preprocess_samples/objects/"
spaceranger_integrated_loc <- paste(objects_loc, "spaceranger.20230816.integrated.rds", sep = "")
spaceranger_metadatad_loc <- paste(objects_loc, "spaceranger.20230816.integrated.arteries.rds", sep = "")

# location of the metadata
metadata_loc <- '/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/metadata/'
tissue_state_loc <- paste(metadata_loc, 'epifat_spot_states.tsv', sep = '')

###############################
# update object
###############################

# read object
spaceranger_integrated <- readRDS(spaceranger_integrated_loc)

# read metadata
tissues_state <- read.table(tissue_state_loc, sep = '\t', header = T)

# add the data
spaceranger_integrated <- add_metadata(seurat_object = spaceranger_integrated, metadata_table = tissues_state)

# save the result
saveRDS(spaceranger_integrated, spaceranger_metadatad_loc)
