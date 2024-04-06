#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: epifat_add_tangram_assignments.R
# Function: add metadata for each section
############################################################################################################################

###############################
# loading libraries
###############################

library(Seurat) #"should be Seurat version 5!"

###############################
# functions
###############################

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
spaceranger_metadatad_loc <- paste(objects_loc, "spaceranger.20230816.integrated.arteries.rds", sep = "")
spaceranger_tangramed_loc <- paste(objects_loc, "spaceranger.20230816.integrated.arteries.fulldecon.rds", sep = "")

# location of tangram output
tangram_tsv_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/deconvolution/tangram/results/'

###############################
# update object
###############################

# read object
spaceranger_integrated <- readRDS(spaceranger_metadatad_loc)

# get the tangram results
tangram_a1_assignements_loc <- paste(tangram_tsv_output_loc, 'epifat_tangram_scores_a1.tsv', sep = '')
tangram_a1_assignments <- read.table(tangram_a1_assignements_loc, sep = '\t', header = T, row.names = 1, check.names = F)
# add the slide to the barcode
rownames(tangram_a1_assignments) <- paste('V10A20_016_A1', rownames(tangram_a1_assignments), sep = '_')
# replace the space with an underscore for the cell types, as we did for RCTD as well
colnames(tangram_a1_assignments) <- gsub(' ', '_', colnames(tangram_a1_assignments))
# now for B to D as well
tangram_b1_assignements_loc <- paste(tangram_tsv_output_loc, 'epifat_tangram_scores_b1.tsv', sep = '')
tangram_b1_assignments <- read.table(tangram_b1_assignements_loc, sep = '\t', header = T, row.names = 1, check.names = F)
rownames(tangram_b1_assignments) <- paste('V10A20_016_B1', rownames(tangram_b1_assignments), sep = '_')
colnames(tangram_b1_assignments) <- gsub(' ', '_', colnames(tangram_b1_assignments))
tangram_c1_assignements_loc <- paste(tangram_tsv_output_loc, 'epifat_tangram_scores_c1.tsv', sep = '')
tangram_c1_assignments <- read.table(tangram_c1_assignements_loc, sep = '\t', header = T, row.names = 1, check.names = F)
rownames(tangram_c1_assignments) <- paste('V10A20_016_C1', rownames(tangram_c1_assignments), sep = '_')
colnames(tangram_c1_assignments) <- gsub(' ', '_', colnames(tangram_c1_assignments))
tangram_d1_assignements_loc <- paste(tangram_tsv_output_loc, 'epifat_tangram_scores_d1.tsv', sep = '')
tangram_d1_assignments <- read.table(tangram_d1_assignements_loc, sep = '\t', header = T, row.names = 1, check.names = F)
rownames(tangram_d1_assignments) <- paste('V10A20_016_D1', rownames(tangram_d1_assignments), sep = '_')
colnames(tangram_d1_assignments) <- gsub(' ', '_', colnames(tangram_d1_assignments))
# merge all four together
tangram_assignments <- rbind(tangram_a1_assignments, tangram_b1_assignments, tangram_c1_assignments, tangram_d1_assignments)

# add as an assay
spaceranger_integrated[['Tangram']] <- CreateAssay5Object(data = t(tangram_assignments[rownames(tangram_assignments) %in% colnames(spaceranger_integrated), ]))

# add a prepend to the column names of tangram cell types, so we can also add them to the metadata
colnames(tangram_assignments) <- paste('tangram', colnames(tangram_assignments), sep = '.')

# now add to the metadata
spaceranger_integrated <- AddMetaData(spaceranger_integrated, metadata = tangram_assignments)

# save the result
saveRDS(spaceranger_integrated, spaceranger_tangramed_loc)
