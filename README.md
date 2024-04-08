# spatial FAT
This repository contains the code that was used for Spatial FAT study, using spatial transcriptomics on epicardial fat tissue. More information can be found in the original paper: TODO

## data availability and description

the dataset consists of one slide, with four sections. section B and D contain one piece of tissue each, however section A and C both contain two pieces of tissue each. Section B and D containg healthy tissue, while sections A and C contain donor healthy and disease tissue from the same donor.

Expression data is available here:
https://eqtlgen.org/sc/datasets/

Expression data is available in two flavours:
- raw spaceranger output
- normalized count matrices
- integrated using CCA as Seurat object


### A. raw spaceranger output

To use the raw spaceranger directories, use the following function and call

Given the spaceranger output was download and present in the spaceranger_out folder, they can be loaded into Seurat, using the following code:
```r
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



###############################
# set the locations of the files
###############################
space_ranger_loc <- './spaceranger_out/'
spaceranger_object_list <- get_spaceranger_objects(space_ranger_loc)

```


## data processing

### alignment to reference genome
alignment was done to the b38 version of the human genome using Spaceranger


### initial QC
'*preprocess/epifat_read_batch1.R*' Spaceranger output was read per section, and QCed using the parameters shown in the script. Then saved in a list per section.


### deconvolution of spots
'*deconvolute/rctd/deconvolution/epifat_create_decon_hca_reference.ipynb*' The 'Heart Global' dataset was downloaded from heartcellatlas.org, and the log-normalized counts were reverted to simple integer counts, then saved in slices of 10k cells for memory efficiency purposes, in .mtx format. The gene list, barcode list, and metadata were also experoted to tsv files.

'*deconvolute/rctd/deconvolution/epifat_created_rctd_reference*' The count matrices and metadata from the previous step were used to create a RCTD spacex reference dataset for deconvolution of spatial data.

'*deconvolute/rctd/deconvolution/epifat_decon_batch1.R*' spacex RCTD is used to deconvolute the spots for each section, with the reference of the previous step being used.

'*deconvolute/rctd/deconvolution/epifat_deconvolution_plots.R*' plot the deconvoluted spots

'*deconvolute/rctd/deconvolution/epifat_spot_proportions.R*' plot the the proportion of celltypes across spots

'*deconvolute/rctd/deconvolution/epifat_pcc_comparison.R*' calculate the correlation of marker gene expression and the proportions of respective cell types per spot


### differential expression based on deconvoluted data
'*deconvolute/rctd/cside/epifat_cside_doublet_mode.R*' do differential expression in A1 and C1, between healthy and sick tissue

'*deconvolute/rctd/cside/epifat_cside_pathways.R*' do pathway enrichment on the DE results of cside DE


### merge sections
'*merge/epifat_merge_batch1.R*' Seurat's CCA integration method is used to integrate the deconvoluted list of Seurat section objects. Then dimensional reduction is performed. Finally the deconvolution data present in the metadata, is added as an assay.


### add metata
'*metadata/epifat_add_metadata_batch1.R*' add metadata on tissue status to the Seurat objectß

### differential expression based on spots
'*conventional_de/epifact_conventional_de.R*' perform DE between clusters calculated over the spots

'*conventional_de/epifat_compare_eat_sat*' compare genes expressed in the slides to ones identified in previous studies

'*conventional_de/epifat_de_dotplot*' plot the conventional DE resulting marker genes


### 

