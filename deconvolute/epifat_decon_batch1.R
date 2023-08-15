#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Irene van Bloklad, Roy Oelen
# Name: epifat_decon_batch1.R
# Function: deconvolute batch1
############################################################################################################################

###############################
# loading libraries
###############################

library(Seurat) #"should be Seurat version 5!"
library(spacexr)
library(Matrix)


###############################
# functions
###############################

prepare_spaceranger_objects_for_integration <- function(spaceranger_object_list, update_seurat_object=F) {
  # check each section
  for (section in names(spaceranger_object_list)) {
    # get the object from the list
    section_object <- spaceranger_object_list[[section]]
    # rename sections
    section_object@meta.data[["section"]] <- gsub("-", "_", section)
    section_object@meta.data[["barcode"]] <- colnames(section_object)
    section_object@meta.data[["section_barcode"]] <- paste(section_object@meta.data[["section"]], section_object@meta.data[["barcode"]], sep = "_")
    section_object <- RenameCells(section_object, new.names = section_object@meta.data[["section_barcode"]])
    
    # update Seurat object if asked
    if (update_seurat_object) {
      section_object <- UpdateSeuratObject(section_object)
    }
    
    # put the object back in the list
    spaceranger_object_list[[section]] <- section_object
  }
  return(spaceranger_object_list)
}



get_decon_per_section <- function(sections, spaceranger_output_loc, reference){
  # save the result per section
  decon_per_section <- list()
  # check each section
  for (section in sections){
    
    # count matrix paths
    rctd_st <- paste(spaceranger_output_loc, section, "/outs/filtered_feature_bc_matrix/", sep = "")
    rctd_st_metadata_loc <- paste(rctd_st, 'matrix.mtx.gz', sep = '')
    rctd_st_genes_loc <- paste(rctd_st, 'features.tsv.gz', sep = '')
    rctd_st_barcodes_loc <- paste(rctd_st, 'barcodes.tsv.gz', sep = '')
    
    # location of the coordinates
    rctd_coordsfile_loc <- paste(spaceranger_output_loc, section, "/outs/spatial/", "tissue_positions.csv", sep = '')
    
    # read the matrix
    rctd_matrix <- readMM(rctd_st_metadata_loc)
    # read the genes
    rctd_genes <- read.table(rctd_st_genes_loc, header = F)$V1
    # read the barcodes
    rctd_barcodes <- read.table(rctd_st_barcodes_loc, header = F)$V1
    # read the coordinate file which is the same for all 4 capture areas
    rctd_coordsfile <- read.table(rctd_coordsfile_loc, header = F)
    
    # organizing the coords matrix
    # extract column names from the first row
    col_names <- unlist(strsplit(rctd_coordsfile[1, 1], ","))
    # remove the first row containing column names
    rctd_coordsfile <- rctd_coordsfile[-1, , drop = FALSE]
    # split the data into separate columns
    rctd_coordsfile_split <- strsplit(rctd_coordsfile$V1, ",", fixed = TRUE)
    rctd_coordsfile_matrix <- do.call(rbind, rctd_coordsfile_split)
    # create a data frame with the split data and set column names
    rctd_coordsfile_df <- as.data.frame(rctd_coordsfile_matrix)
    colnames(rctd_coordsfile_df) <- col_names
    
    # generating the counts matrix (rows = genes, cols = barcodes)
    # set the barcodes and genes
    colnames(rctd_matrix) <- rctd_barcodes
    rownames(rctd_matrix) <- rctd_genes
    
    # generating the coords matrix (rows = barcode, cols = x and y)
    # array_col is x-coordinate and array_row is y-coordinate
    coords <- data.frame(
      x = as.numeric(rctd_coordsfile_df$array_col),
      y = as.numeric(rctd_coordsfile_df$array_row)
    )
    rownames(coords) <- rctd_coordsfile_df$barcode
    
    # generating nUMI which is the total counts per spot
    nUMI <- colSums(rctd_matrix)
    
    ### Create SpatialRNA object
    puck <- SpatialRNA(coords, rctd_matrix, nUMI)
    
    # pixels/spots to be used (a list of barcode names). 
    barcodes <- colnames(puck@counts)
    
    # creating RCTD Object
    myRCTD <- create.RCTD(puck, reference, max_cores = 1)
    
    # running RCTD
    myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
    
    # RCTD results
    results <- myRCTD@results
    
    # put result in list
    decon_per_section[[section]] <- results
  }
  return(decon_per_section)
}


add_decons_to_sections <- function(spaceranger_object_list, decon_per_section){
  # we can only do the sections we have seurat and decon info for
  sections <- intersect(names(spaceranger_object_list), names(decon_per_section))
  # check each section
  for (section in sections) {
    # grab the result
    results <- decon_per_section[[section]]
    # normalize the cell type proportions to sum to 1.
    norm_weights = normalize_weights(results$weights) 
    
    # get the object from the list
    section_object <- spaceranger_object_list[[section]]
    
    # add the deconvolution result to the object
    # extract the weights from RCTD
    norm_weights <- t(norm_weights)
    # add the slide to the colnames
    colnames(norm_weights) <- paste(gsub("-", "_", section), colnames(norm_weights), sep = "_")
    # get the cells that are present in the seurat objects
    cells_seurat <- colnames(section_object)
    # get the cells that are present in the decon object
    cells_decon <- colnames(norm_weights)
    # get the cells that are present in both
    cells_both <- intersect(cells_seurat, cells_decon)
    # warn if there are cells missing in either
    cells_only_seurat <- setdiff(cells_seurat, cells_both)
    cells_only_decon <- setdiff(cells_decon, cells_both)
    if (length(cells_only_decon) > 0) {
      message(paste("dropping", length(cells_only_decon), "cells from decon that are not present in Seurat"))
      print(head(cells_only_decon))
    }
    if (length(cells_only_seurat) > 0) {
      message(paste("dropping", length(cells_only_seurat), "cells from Seurat that are not present in decon"))
      print(head(cells_only_seurat))
    }
    # keep only the barcodes we have in both
    norm_weights[, cells_both]
    section_object <- section_object[, cells_both]
    
    # create new assay and object
    # rctd_assay <- CreateAssay5Object(data = norm_weights)
    
    # and add it to the original object
    # section_object[["RCTD"]] <- rctd_assay
    
    # add as metadata
    section_object <- AddMetaData(section_object, metadata = t(norm_weights))

    # put the object back in the list
    spaceranger_object_list[[section]] <- section_object
  }
  return(spaceranger_object_list)
}

###############################
# set some parameters
###############################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 2000 * 1000 * 1024^2)

# set seed
set.seed(7777)


###############################
# set the locations of the files
###############################
space_ranger_loc <- '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/epicardial_fat/processed/alignment/spaceranger_out/'

# location to store the resulting files
objects_loc <- "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/seurat_preprocess_samples/objects/"
spaceranger_object_loc <- paste(objects_loc, "spaceranger.20230811.rds", sep = "")
spaceranger_deconvolute_loc <- paste(objects_loc, "spaceranger.20230811.deconvoluted.rds", sep = "")

# read the object
spaceranger_object_list <- readRDS(spaceranger_object_loc)

# prepare for SCT integration
spaceranger_object_list <- prepare_spaceranger_objects_for_integration(spaceranger_object_list, update_seurat_object = T)

###############################
# 5. Deconvolution - per section using raw counts
###############################

# directory of the reference dataset which Roy prepared
reference <- readRDS("/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/epicardial_fat/ongoing/rtcd/references/hca/HCA_RTCD_ref.rds")

# do deconvolution
decon_per_section <- get_decon_per_section(names(spaceranger_object_list), space_ranger_loc, reference)

# add to the objects
spaceranger_object_list <- add_decons_to_sections(spaceranger_object_list, decon_per_section)

# save the result
saveRDS(spaceranger_object_list, spaceranger_deconvolute_loc)
