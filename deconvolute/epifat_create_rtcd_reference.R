#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: epifat_create_rtcd_reference.R
# Function: create the RTCD reference dataset using the sliced count matrices
############################################################################################################################


####################
# libraries        #
####################

library(Matrix) # to read the matrix slices
library(spacexr) # to create the reference

####################
# Functions        #
####################

#' read a directory containing .mtx files that were sliced into one sparse matrix
#' 
#' @param slices_loc the directory containing the .mtx files
#' @param slices_regex the regex for the .mtx files in the directory
#' @param verbose whether or not to to mention which slice we are at
#' @returns a sparse matrix of all the combined slices
#' merged_slices <- read_sliced_matrices(hca_loc)
read_sliced_matrices <- function(slices_loc, slices_regex='*.mtx', slice_sort='(\\d+)', verbose=F) {
  # make a list to store the slices
  slices_list <- list()
  # do the listing of files
  slices <- list.files(slices_loc, pattern = slices_regex, full.names = F)
  # sort the slices if we need to
  if (!is.null(slice_sort)) {
    # get the number regex returns
    slice_regex_returns <- regexpr(slice_sort, slices)
    # store the slice numbers
    slice_starts <- rep(NA, times = length(slice_regex_returns))
    # check each slice regex
    for (regex_i in 1:length(slice_regex_returns)) {
      # extract the start of the regex match
      regex_start <- slice_regex_returns[[regex_i]][[1]]
      # and the length of the regex
      regex_length <- attr(x = slice_regex_returns, which = 'match.length')[regex_i]
      # extract the start from the filename
      slice_start_string <- substr(slices[regex_i], regex_start, regex_start + regex_length - 1)
      # add to the slice starts
      slice_starts[regex_i] <- as.numeric(slice_start_string)
    }
    # now order the slices based on these slice starts
    slices <- slices[order(slice_starts)] 
  }
  # check each slice
  for (slice_loc in slices) {
    # create the full path
    slice_path_full <- paste(slices_loc, slice_loc, sep = '')
    # read the slice
    slice <- readMM(slice_path_full)
    # transpose it, so that it is genes on the rows, and cells on the column
    slice <- t(slice)
    # store in the list
    slices_list[[slice_loc]] <- slice
  }
  # merge the slices
  merged_slices <- do.call('cbind', slices_list)
  return(merged_slices)
}

#' replace a given vector of celltypes illegal characters with posix safe ones
#' 
#' @param cell_types the vector of cell types to make safe
#' @returns a vector of the same cell types, but posix safe
#' cell_types <- make_celltypes_safe(hca_metadata$cell_type)
make_celltypes_safe <- function(cell_types){
  # get a safe file name
  cell_type_safes <- gsub(' |/', '_', cell_types)
  cell_type_safes <- gsub('-', '_negative', cell_type_safes)
  cell_type_safes <- gsub('\\+', '_positive', cell_type_safes)
  cell_type_safes <- gsub('\\)', '', cell_type_safes)
  cell_type_safes <- gsub('\\(', '', cell_type_safes)
  return(cell_type_safes)
}


####################
# Main Code        #
####################

# location of the data
hca_loc <- '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/epicardial_fat/ongoing/rtcd/references/hca/'
hca_metadata_loc <- paste(hca_loc, 'hca_metadata.tsv', sep = '')
hca_genes_loc <- paste(hca_loc, 'hca_genes.tsv', sep = '')
hca_barcodes_loc <- paste(hca_loc, 'hca_barcodes.tsv', sep = '')

# rtcd reference location
rctd_ref_loc <- paste(hca_loc, 'HCA_RTCD_ref.rds', sep = '')

# read the metadata
hca_metadata <- read.table(hca_metadata_loc, sep = '\t', header = T, row.names = 1)
# read the genes
hca_genes <- read.table(hca_genes_loc, header = F)$V1
# read the barcodes
hca_barcodes <- read.table(hca_barcodes_loc, header = F)$V1

# read all the slices
merged_slices <- read_sliced_matrices(hca_loc)
# set the barcodes and genes
rownames(merged_slices) <- hca_genes
colnames(merged_slices) <- hca_barcodes

# now do the RTCD reference preparation
cell_types <- make_celltypes_safe(hca_metadata$cell_type) # get the cell types, but without illegal characters
names(cell_types) <- hca_barcodes # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- hca_metadata$total_counts # get the total counts from the metadata
names(nUMI) <- hca_barcodes # create nUMI named list
reference <- Reference(merged_slices, cell_types, nUMI) # create the object
saveRDS(reference, rctd_ref_loc) # save result
