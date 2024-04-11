# spatial FAT
This repository contains the code that was used for Spatial FAT study, using spatial transcriptomics on epicardial fat tissue. More information can be found in the original paper

## data availability and description

The dataset consists of one slide, with four sections. Sections B and D contain one piece of tissue each, however sections A and C both contain two pieces of tissue each. Section B and D containg healthy tissue, while sections A and C contain donor healthy and disease tissue from the same donor.

Expression data and more information is available here:
https://downloads.molgeniscloud.org/downloads/epifat/

Expression data is available in two flavours:
- raw spaceranger output (https://downloads.molgeniscloud.org/downloads/epifat/epifat_spaceranger_outs.zip)
- integrated using CCA as Seurat object (https://downloads.molgeniscloud.org/downloads/epifat/epifat_seurat_integrated_decon_metadata.rds)


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

# read objects
space_ranger_loc <- './spaceranger_out/'
spaceranger_object_list <- get_spaceranger_objects(space_ranger_loc)

```

Or in python like this
```python
def read_slices(slices_loc, slices, counts_file='filtered_feature_bc_matrix.h5', do_norm=True, do_dimruc=True):
    """read ST slices, do QC, and put them in a list
        
        Parameters
        ----------
        slices_loc : str
            the location of the folders containing the slices
        slices : list
            a list containing the slices (directory names)
        counts_file : str
            the name of the count expression matrix
        do_norm : bool
            run normalization
        do_dimruc : bool
            run dimensional reduction
        
        Returns
        -------
        result
           a dictionary of AnnData objects
        """
    # create a dictionary to store the slices
    slices_dict = {}
    # read each slice
    for slice in slices:
        # paste together the path
        full_visium_path = ''.join([slices_loc, '/', slice, '/outs/'])
        # read the file
        adata = sc.read_visium(path = full_visium_path,
                             count_file = counts_file)
        # make gene names unique
        adata.var_names_make_unique()
        # do some standard preprocessing
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

        sc.pp.filter_cells(adata, min_counts=200)
        sc.pp.filter_genes(adata, min_cells=3)
        
        # normalize if requested
        if do_norm:
            # as well as the normalization
            sc.pp.normalize_total(adata, inplace=True)
            sc.pp.log1p(adata)
            # and calculation of highly variable genes
            sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
            
            # and dimensional reduction if requested
            if do_dimruc:
                # calculate principal components as well
                sc.pp.pca(adata)
                # do nearest neighbours
                sc.pp.neighbors(adata)
                # do 2d UMAP dim reduction
                sc.tl.umap(adata)
                # and Leiden clustering
                sc.tl.leiden(adata, key_added="clusters")
        
        # add the result to the dictionary
        slices_dict[slice] = adata
    return slices_dict
    
# read each slices
slice_objects = read_slices('./spaceranger_out/', ['V10A20-016_A1', 'V10A20-016_B1', 'V10A20-016_C1', 'V10A20-016_D1'])

```


### B. Seurat object

The Seurat object can simply be read using Seurat in R

```r

epifat <- readRDS('epifat_seurat_integrated_decon_metadata.rds')

```


This Rstudio image with included libraries was used to process the data: https://github.com/royoelen/single-cell-container-server

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
