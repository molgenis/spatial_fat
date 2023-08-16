# spatial_fat
spatial FAT study, using spatial transcriptomics on epicardial fat tissue

## dataset
the dataset consists of one slide, with four sections. section B and D contain one piece of tissue each, however section A and C both contain two pieces of tissue each. Section B and D containg healthy tissue, while sections A and C contain donor healthy and disease tissue from the same donor.

## alignment to reference genome
alignment was done to the b38 version of the human genome using Spaceranger

## initial QC
'*preprocess/epifat_read_batch1.R*' Spaceranger output was read per section, and QCed using the parameters shown in the script. Then saved in a list per section.

## deconvolution of spots
'*deconvolute/epifat_create_decon_hca_reference.ipynb*' The 'Heart Global' dataset was downloaded from heartcellatlas.org, and the log-normalized counts were reverted to simple integer counts, then saved in slices of 10k cells for memory efficiency purposes, in .mtx format. The gene list, barcode list, and metadata were also experoted to tsv files.

'*deconvolute/epifat_created_rctd_reference*' The count matrices and metadata from the previous step were used to create a RCTD spacex reference dataset for deconvolution of spatial data.

'*deconvolute/epifat_decon_batch1.R*' spacex RCTD is used to deconvolute the spots for each section, with the reference of the previous step being used.

## merge sections
'*merge/epifat_merge_batch1*' Seurat's CCA integration method is used to integrate the deconvoluted list of Seurat section objects. Then dimensional reduction is performed. Finally the deconvolution data present in the metadata, is added as an assay.