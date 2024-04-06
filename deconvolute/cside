############################################################################################################################
# Authors: Irene van Blokland, Roy Oelen, Aida
# Name: epifat_cside.R
# Function: Cell type-specific DE, accounting for changes in cell type proportions and cell mixtures
############################################################################################################################

###############################
# loading libraries
###############################
library(Seurat)
library(spacexr)
library(Matrix)
library(doParallel)
library(ggplot2)
library(biomaRt)

###############################
# getting correct Seurat version
###############################
# RCTD should be performed in full mode
# RCTD results are stored in @results$weights where each entry represents the estimated proportion of each cell type on each spot

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 2000 * 1000 * 1024^2)

# set seed
set.seed(7777)

###############################
# RCTD functions
###############################

get_decon_per_section_doublet <- function(sections, spaceranger_output_loc, reference){
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
    myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
    
    # put result in list
    decon_per_section[[section]] <- myRCTD
  }
  return(decon_per_section)
}

###############################
# Do RCTD 
###############################
space_ranger_loc <- "/groups/umcg-franke-scrna/tmp02/projects/epifat/processed/alignment/spaceranger_out/"

# location to store the resulting files
objects_loc <- "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/seurat_preprocess_samples/objects/"
spaceranger_object_loc <- paste(objects_loc, "spaceranger.20230811.rds", sep = "")
spaceranger_deconvolute_loc <- paste(objects_loc, "spaceranger.20230811.deconvoluted.rds", sep = "")
spaceranger_rctd_loc_d <- paste(objects_loc, "spaceranger.20230912.rctd.doublet.rds", sep = "")
spaceranger_status_loc <- paste(objects_loc, "spaceranger.20230816.integrated.arteries.rds", sep = "")

### perform RCTD just to generate the object
# read the object
spaceranger_object_list <- readRDS(spaceranger_object_loc)
spaceranger_object <- readRDS(spaceranger_status_loc)
# directory of the reference dataset which Roy prepared
reference <- readRDS("/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/epicardial_fat/ongoing/rtcd/references/hca/HCA_RTCD_ref.rds")
# do deconvolution
decon_per_section_all_d <- get_decon_per_section_doublet(names(spaceranger_object_list), space_ranger_loc, reference)
# save the RCTD object
saveRDS(decon_per_section_all_d, file = "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/seurat_preprocess_samples/objects/spaceranger.20230912.rctd.doublet.rds")

###############################
# Prepare for C-SIDE, patient A1
###############################
# read the deconvoluted raw object
spaceranger_rctd_list_d <- readRDS(spaceranger_rctd_loc_d)
# select section V10A20-016_A1
myRCTD_section_A1 <- spaceranger_rctd_list_d[["V10A20-016_A1"]]
# make a subset with barcodes as rownames and celltype proportions as columns
barcodes <- colnames(myRCTD_section_A1@spatialRNA@counts)
weights <- myRCTD_section_A1@results$weights
norm_weights <- normalize_weights(weights)
cell_types <- c('Adipocyte','Atrial_Cardiomyocyte', "Endothelial_cell", "Fibroblast", "Lymphatic_Endothelial_cell", "Lymphoid", "Mast_cell", "Mesothelial_cell", "Mural_cell", "Myeloid", "Neural_cell", "Ventricular_Cardiomyocyte")
print(head(norm_weights[,cell_types]))

######### Trying this method: https://raw.githack.com/dmcable/spacexr/master/vignettes/visium_full_regions.html
## choose regions of control and disease, fetching the data from the integrated object
# load barcodes from other object
status <- readRDS(spaceranger_status_loc)
# select only A1
status_A1 <- subset(status, section == "V10A20_016_A1")
barcodes <- status_A1$barcode
status <- status_A1$status
# status_vector <- paste(barcodes, status, sep = "_")
barcodes_S <- unname(status_A1$barcode[status_A1$status == "S"])
barcodes_C <- unname(status_A1$barcode[status_A1$status == "H"])
# making the region list
region_list = list("1" = barcodes_S, "2" = barcodes_C)
# removing the names of the list
explanatory.variable <- unname(region_list)

####### trying this method: https://github.com/dmcable/spacexr/blob/master/vignettes/CSIDE_two_regions.Rmd
# for this method, the explanatory var is an integer with barcodes as names and 0 or 1 as value. 
# so Healthy = 0, sick = 1. 

#### sorting the spots in the rctd object to the healthy/disease spots
# subsetting the rctd object
coords <- myRCTD_section_A1@spatialRNA@coords
counts <- myRCTD_section_A1@spatialRNA@counts
umi <- myRCTD_section_A1@spatialRNA@nUMI

#### making a list of barcodes to remove from the RCTD object
# selecting spots that are in access in the integrated object compared to the RCTD
bar_rctd <- colnames(myRCTD_section_A1@spatialRNA@counts) # length 1219
bar_int <- status_A1@meta.data$barcode # length 1112
# identifying the differing barcodes
bar_dif <- setdiff(bar_rctd, bar_int) 
# removing the identified barcodes from the rctd object
bar_rctd_t <- setdiff(bar_rctd, bar_dif) # length 1112
# matching the order of the integrated object to the rctd object
bc_counts <- colnames(myRCTD_section_A1@spatialRNA@counts)
bc_order <- bc_counts
# setting control barcodes to 0, disease to 1
H_vec <- rep(0, length(barcodes_C))
names(H_vec) <- barcodes_C
S_vec <- rep(1, length(barcodes_S))
names(S_vec) <- barcodes_S
# combining the two
exp_var <- c(H_vec, S_vec)
# fetching barcodes
exp_var.names <- names(exp_var)
# sorting the order
exp_var.t <- exp_var[match(bc_order, exp_var.names)]
# checking the order, making sure it is the same
bar_rctd.t <- bar_rctd[match(bar_rctd_t, bar_rctd)]
#identical(exp_var.t, bc_vec.names)
# getting barcodes from all three parts 
bc_counts <- colnames(myRCTD_section_A1@spatialRNA@counts)
bc_coords <- rownames(myRCTD_section_A1@spatialRNA@coords)
bc_numi <- names(myRCTD_section_A1@spatialRNA@nUMI)
# checking whether barcodes are identical
identical(bc_coords, bc_counts) #TRUE
identical(bc_coords, bc_numi) #TRUE
identical(bc_counts, bc_numi) #TRUE
identical(bc_counts, bc_coords) #TRUE
identical(bc_counts, bc_numi) #TRUE
identical(bc_coords, bc_numi) #TRUE
# filtering the barcodes of the spatial object in the right order
myRCTD_section_A1@spatialRNA@coords <- myRCTD_section_A1@spatialRNA@coords[rownames(myRCTD_section_A1@spatialRNA@coords)%in%bar_int,]
myRCTD_section_A1@spatialRNA@nUMI <- myRCTD_section_A1@spatialRNA@nUMI[names(myRCTD_section_A1@spatialRNA@nUMI)%in%bar_int]
myRCTD_section_A1@spatialRNA@counts <- myRCTD_section_A1@spatialRNA@counts[,bar_int]

# sorting the explanatory variable so the order matches the RCTD object
sorting_index <- order(bc_coords)
exp_var_reorder <- explanatory.variable[bc_coords]
explanatory.variable <- c(H_vec, S_vec)
explanatory.variable <- explanatory.variable[match(bc_counts, names(explanatory.variable))]

###############################
# Choose regions/generate explanatory variable 
###############################
# Each covariate should be a numeric vector,  representing the value of the 
# covariate at each pixel. Please standardize each covariate between 0 and 1, and the names of each 
# covariate variable should match the names, or barcodes, of the spatial transcriptomics pixels.

# lets try using clusters as regions
#barcodes_1 <- myRCTD_section$barcode[myRCTD_section$seurat_clusters == "1"]
#barcodes_2 <- myRCTD_section$barcode[myRCTD_section$seurat_clusters == "2"]
#barcodes_3 <- myRCTD_section$barcode[myRCTD_section$seurat_clusters == "3"]
#barcodes_4 <- myRCTD_section$barcode[myRCTD_section$seurat_clusters == "4"]
#barcodes_5 <- myRCTD_section$barcode[myRCTD_section$seurat_clusters == "5"]
#barcodes_6 <- myRCTD_section$barcode[myRCTD_section$seurat_clusters == "6"]
#barcodes_7 <- myRCTD_section$barcode[myRCTD_section$seurat_clusters == "7"]
#barcodes_8 <- myRCTD_section$barcode[myRCTD_section$seurat_clusters == "8"]
#barcodes_9 <- myRCTD_section$barcode[myRCTD_section$seurat_clusters == "9"]
#barcodes_10 <- myRCTD_section$barcode[myRCTD_section$seurat_clusters == "10"]
#barcodes_11 <- myRCTD_section$barcode[myRCTD_section$seurat_clusters == "11"]
#barcodes_12 <- myRCTD_section$barcode[myRCTD_section$seurat_clusters == "12"]
#barcodes_13 <- myRCTD_section$barcode[myRCTD_section$seurat_clusters == "13"]
#barcodes_14 <- myRCTD_section$barcode[myRCTD_section$seurat_clusters == "14"]

#region_list = list(barcodes_1, barcodes_2, barcodes_3, barcodes_4, barcodes_5, barcodes_6,
#barcodes_7, barcodes_8, barcodes_9, barcodes_10, barcodes_11, barcodes_12,
#barcodes_13, barcodes_14)

#explanatory.variable <- c(myRCTD_section$barcode)
#names(explanatory.variable) <- myRCTD_section$status
#names(explanatory.variable) <- ifelse(myRCTD_section$status == "H", 0, 1)
# select barcodes to match state
# barcodes_S <- myRCTD_section$barcode[myRCTD_section$status == "S"]
# barcodes_C <- myRCTD_section$barcode[myRCTD_section$status == "H"]
# create list of regions
# Create a vector of all unique barcodes
# all_barcodes <- unique(c(barcodes_S, barcodes_C))
# Initialize an atomic vector to store the results
# barcode_vector <- data.frame(
#  barcode = all_barcodes,
#  H = as.integer(all_barcodes %in% barcodes_C),  # 1 if in barcodes_C, 0 otherwise
#  S = as.integer(all_barcodes %in% barcodes_S)   # 1 if in barcodes_S, 0 otherwise
#)
# set barcode vector to explanatory variable
# explanatory.variable <- barcode_vector

###############################
# Running C-SIDE, patient A1
###############################

# run CSIDE
myRCTD_section_A1 <- run.CSIDE.single(myRCTD_section_A1, explanatory.variable, gene_threshold = .00005, cell_type_threshold = 125, fdr = 0.01, weight_threshold = NULL)
# save results
saveRDS(myRCTD_section_A1, "/groups/umcg-franke-scrna/tmp01/projects/epifat/ongoing/seurat_preprocess_samples/objects/myRCTD_section_A1.230912.cside.rds")

# load cside
myRCTD_section_A1 <- readRDS("/groups/umcg-franke-scrna/tmp01/projects/epifat/ongoing/seurat_preprocess_samples/objects/myRCTD_section_A1.230912.cside.rds")

###############################
# Explore CSIDE results, patient A1
###############################
# significant genes are stored in myRCTD@de_results$sig_gene_list
results_de_A1 <- myRCTD_section_A1@de_results$sig_gene_list$Adipocyte
results_de_A1 <- myRCTD_section_A1@de_results$sig_gene_list$Atrial_Cardiomyocyte
results_de_A1 <- myRCTD_section_A1@de_results$sig_gene_list$Endothelial_cell
results_de_A1 <- myRCTD_section_A1@de_results$sig_gene_list$Fibroblast
results_de_A1 <- myRCTD_section_A1@de_results$sig_gene_list$Lymphatic_Endothelial_cell
results_de_A1 <- myRCTD_section_A1@de_results$sig_gene_list$Lymphoid
results_de_A1 <- myRCTD_section_A1@de_results$sig_gene_list$Mast_cell
results_de_A1 <- myRCTD_section_A1@de_results$sig_gene_list$Mesothelial_cell
results_de_A1 <- myRCTD_section_A1@de_results$sig_gene_list$Mural_cell
results_de_A1 <- myRCTD_section_A1@de_results$sig_gene_list$Myeloid
results_de_A1 <- myRCTD_section_A1@de_results$sig_gene_list$Neural_cell
results_de_A1 <- myRCTD_section_A1@de_results$sig_gene_list$Ventricular_Cardiomyocyte

# make a list of all DE genes for section A1
results_de_A1_adipo <- myRCTD_section_A1@de_results$all_gene_list$Adipocyte
results_de_A1_endo <- myRCTD_section_A1@de_results$all_gene_list$Endothelial_cell
# save DE lists for section A1
write.csv(results_de_A1_adipo, "/groups/umcg-franke-scrna/tmp02/users/umcg-ivanblokland/cside_doublet_A1_adipo_DE.csv")
write.csv(results_de_A1_endo, "/groups/umcg-franke-scrna/tmp02/users/umcg-ivanblokland/cside_doublet_A1_endo_DE.csv")

# exploring results
dataframe_A1 <- data.frame(results_de_A1)
cell_type = "Endothelial_cell"
sig_gene <- "ENSG00000184557"
# plotting the spatial visualisation of a DE gene
plot_gene_two_regions(myRCTD_section_A1, sig_gene, cell_type, min_UMI = 10)

###############################
# Prepare C-SIDE, patient C1
###############################

# select section V10A20-016_C1
myRCTD_section_C1 <- spaceranger_rctd_list_d[["V10A20-016_C1"]]
# make a subset with barcodes as rownames and celltype proportions as columns
barcodes <- colnames(myRCTD_section_C1@spatialRNA@counts)
weights <- myRCTD_section_C1@results$weights
norm_weights <- normalize_weights(weights)
cell_types <- c('Adipocyte','Atrial_Cardiomyocyte', "Endothelial_cell", "Fibroblast", "Lymphatic_Endothelial_cell", "Lymphoid", "Mast_cell", "Mesothelial_cell", "Mural_cell", "Myeloid", "Neural_cell", "Ventricular_Cardiomyocyte")
print(head(norm_weights[,cell_types]))

######### Trying this method: https://raw.githack.com/dmcable/spacexr/master/vignettes/visium_full_regions.html
## choose regions of control and disease, fetching the data from the integrated object
# load barcodes from other object
status <- readRDS(spaceranger_status_loc)
# select only A1
status_C1 <- subset(status, section == "V10A20_016_C1")
barcodes <- status_C1$barcode
status <- status_C1$status
# status_vector <- paste(barcodes, status, sep = "_")
barcodes_S <- unname(status_C1$barcode[status_C1$status == "S"])
barcodes_C <- unname(status_C1$barcode[status_C1$status == "H"])
# making the region list
region_list = list("1" = barcodes_S, "2" = barcodes_C)
# removing the names of the list
explanatory.variable <- unname(region_list)

####### trying this method: https://github.com/dmcable/spacexr/blob/master/vignettes/CSIDE_two_regions.Rmd
# for this method, the explanatory var is an integer with barcodes as names and 0 or 1 as value. 
# so Healthy = 0, sick = 1. 

#### sorting the spots in the rctd object to the healthy/disease spots
# subsetting the rctd object
coords <- myRCTD_section_C1@spatialRNA@coords
counts <- myRCTD_section_C1@spatialRNA@counts
umi <- myRCTD_section_C1@spatialRNA@nUMI

#### making a list of barcodes to remove from the RCTD object
# selecting spots that are in access in the integrated object compared to the RCTD
bar_rctd <- colnames(myRCTD_section_C1@spatialRNA@counts) # length 1219
bar_int <- status_C1@meta.data$barcode # length 1112
# identifying the differing barcodes
bar_dif <- setdiff(bar_rctd, bar_int) 
# removing the identified barcodes from the rctd object
bar_rctd_t <- setdiff(bar_rctd, bar_dif) # length 1112
# matching the order of the integrated object to the rctd object
bc_counts <- colnames(myRCTD_section_C1@spatialRNA@counts)
bc_order <- bc_counts
# setting control barcodes to 0, disease to 1
H_vec <- rep(0, length(barcodes_C))
names(H_vec) <- barcodes_C
S_vec <- rep(1, length(barcodes_S))
names(S_vec) <- barcodes_S
# combining the two
exp_var <- c(H_vec, S_vec)
# fetching barcodes
exp_var.names <- names(exp_var)
# sorting the order
exp_var.t <- exp_var[match(bc_order, exp_var.names)]
# checking the order, making sure it is the same
bar_rctd.t <- bar_rctd[match(bar_rctd_t, bar_rctd)]
#identical(exp_var.t, bc_vec.names)
# getting barcodes from all three parts 
bc_counts <- colnames(myRCTD_section_C1@spatialRNA@counts)
bc_coords <- rownames(myRCTD_section_C1@spatialRNA@coords)
bc_numi <- names(myRCTD_section_C1@spatialRNA@nUMI)
# checking whether barcodes are identical
identical(bc_coords, bc_counts) #TRUE
identical(bc_coords, bc_numi) #TRUE
identical(bc_counts, bc_numi) #TRUE
identical(bc_counts, bc_coords) #TRUE
identical(bc_counts, bc_numi) #TRUE
identical(bc_coords, bc_numi) #TRUE
# filtering the barcodes of the spatial object in the right order
myRCTD_section_C1@spatialRNA@coords <- myRCTD_section_C1@spatialRNA@coords[rownames(myRCTD_section_C1@spatialRNA@coords)%in%bar_int,]
myRCTD_section_C1@spatialRNA@nUMI <- myRCTD_section_C1@spatialRNA@nUMI[names(myRCTD_section_C1@spatialRNA@nUMI)%in%bar_int]
myRCTD_section_C1@spatialRNA@counts <- myRCTD_section_C1@spatialRNA@counts[,bar_int]

# sorting the explanatory variable so the order matches the RCTD object
sorting_index <- order(bc_coords)
exp_var_reorder <- explanatory.variable[bc_coords]
explanatory.variable <- c(H_vec, S_vec)
explanatory.variable <- explanatory.variable[match(bc_counts, names(explanatory.variable))]

###############################
# Running C-SIDE, patient C1
###############################
# run CSIDE
myRCTD_section_C1 <- run.CSIDE.single(myRCTD_section_C1, explanatory.variable, gene_threshold = .00005, cell_type_threshold = 125, fdr = 0.01, weight_threshold = NULL)
# save results
saveRDS(myRCTD_section_C1, "/groups/umcg-franke-scrna/tmp01/projects/epifat/ongoing/seurat_preprocess_samples/objects/myRCTD_section_C1.230913.cside.rds")

# load cside
myRCTD_section_C1 <- readRDS("/groups/umcg-franke-scrna/tmp01/projects/epifat/ongoing/seurat_preprocess_samples/objects/myRCTD_section_C1.230913.cside.rds")

###############################
# Explore CSIDE results, patients A1 and C1
###############################
# significant genes are stored in myRCTD@de_results$sig_gene_list
results_de_A1_adipocyte <- myRCTD_section_A1@de_results$sig_gene_list$Adipocyte
results_de_A1_acm <- myRCTD_section_A1@de_results$sig_gene_list$Atrial_Cardiomyocyte # No DE genes
results_de_A1_endo <- myRCTD_section_A1@de_results$sig_gene_list$Endothelial_cell
results_de_A1_fib <- myRCTD_section_A1@de_results$sig_gene_list$Fibroblast # No DE genes
results_de_A1_lymphendo <- myRCTD_section_A1@de_results$sig_gene_list$Lymphatic_Endothelial_cell # No DE genes
results_de_A1_lymph <- myRCTD_section_A1@de_results$sig_gene_list$Lymphoid # No DE genes
results_de_A1_mast <- myRCTD_section_A1@de_results$sig_gene_list$Mast_cell # No DE genes
results_de_A1_meso <- myRCTD_section_A1@de_results$sig_gene_list$Mesothelial_cell # No DE genes
results_de_A1_mural <- myRCTD_section_A1@de_results$sig_gene_list$Mural_cell # No DE genes
results_de_A1_mye <- myRCTD_section_A1@de_results$sig_gene_list$Myeloid # No DE genes
results_de_A1_neu <- myRCTD_section_A1@de_results$sig_gene_list$Neural_cell # No DE genes
results_de_A1_vcm <- myRCTD_section_A1@de_results$sig_gene_list$Ventricular_Cardiomyocyte # No DE genes

results_de_C1_adipocyte <- myRCTD_section_C1@de_results$sig_gene_list$Adipocyte
results_de_C1_acm <- myRCTD_section_C1@de_results$sig_gene_list$Atrial_Cardiomyocyte # No DE genes
results_de_C1_endo <- myRCTD_section_C1@de_results$sig_gene_list$Endothelial_cell
results_de_C1_fib <- myRCTD_section_C1@de_results$sig_gene_list$Fibroblast # No DE genes
results_de_C1_lymphendo <- myRCTD_section_C1@de_results$sig_gene_list$Lymphatic_Endothelial_cell # No DE genes
results_de_C1_lymph <- myRCTD_section_C1@de_results$sig_gene_list$Lymphoid # No DE genes
results_de_C1_mast <- myRCTD_section_C1@de_results$sig_gene_list$Mast_cell # No DE genes
results_de_C1_meso <- myRCTD_section_C1@de_results$sig_gene_list$Mesothelial_cell # No DE genes
results_de_C1_mural <- myRCTD_section_C1@de_results$sig_gene_list$Mural_cell # No DE genes
results_de_C1_mye <- myRCTD_section_C1@de_results$sig_gene_list$Myeloid # No DE genes
results_de_C1_neu <- myRCTD_section_C1@de_results$sig_gene_list$Neural_cell # No DE genes
results_de_C1_vcm <- myRCTD_section_C1@de_results$sig_gene_list$Ventricular_Cardiomyocyte # No DE genes

# make a list of all DE genes for section C1
results_de_C1_adipo <- myRCTD_section_C1@de_results$all_gene_list$Adipocyte
results_de_C1_endo <- myRCTD_section_C1@de_results$all_gene_list$Endothelial_cell

# save DE lists for section A1
write.csv(results_de_C1_adipo, "/groups/umcg-franke-scrna/tmp01/users/umcg-ivanblokland/cside_doublet_C1_adipo_DE.csv")
write.csv(results_de_C1_endo, "/groups/umcg-franke-scrna/tmp01/users/umcg-ivanblokland/cside_doublet_C1_endo_DE.csv")

# make a scatterplot
A1_adipo <- data.frame(Gene = rownames(myRCTD_section_A1@de_results$all_gene_list$Adipocyte), log_fc = myRCTD_section_A1@de_results$all_gene_list$Adipocyte$log_fc)
C1_adipo <- data.frame(Gene = rownames(myRCTD_section_C1@de_results$all_gene_list$Adipocyte), log_fc = myRCTD_section_C1@de_results$all_gene_list$Adipocyte$log_fc)
A1_endo <- data.frame(Gene = rownames(myRCTD_section_A1@de_results$all_gene_list$Endothelial_cell), log_fc = myRCTD_section_A1@de_results$all_gene_list$Endothelial_cell$log_fc)
C1_endo <- data.frame(Gene = rownames(myRCTD_section_C1@de_results$all_gene_list$Endothelial_cell), log_fc = myRCTD_section_C1@de_results$all_gene_list$Endothelial_cell$log_fc)

# merge dataframes per celltype
merged_adipo <- merge(A1_adipo, C1_adipo, by = "Gene", suffixes = c("_A1_adipo", "_C1_adipo"), all = TRUE)
merged_endo <- merge(A1_endo, C1_endo, by = "Gene", suffixes = c("_A1_endo", "_C1_endo"), all = TRUE)

# Create a new column for gene specificity
merged_adipo <- merged_adipo %>%
  mutate(
    specificity = case_when(
      is.na(log_fc_C1_adipo) ~ "A1_specific",
      is.na(log_fc_A1_adipo) ~ "C1_specific",
      TRUE ~ "Matching"
    )
  )
merged_endo <- merged_endo %>%
  mutate(
    specificity = case_when(
      is.na(log_fc_C1_endo) ~ "A1_specific",
      is.na(log_fc_A1_endo) ~ "C1_specific",
      TRUE ~ "Matching"
    )
  )
  
# save the merged files
write.csv(merged_adipo, "/groups/umcg-franke-scrna/tmp02/users/umcg-ivanblokland/cside_doublet_adipo_DE.csv")
write.csv(merged_endo, "/groups/umcg-franke-scrna/tmp02/users/umcg-ivanblokland/cside_doublet_endo_DE.csv")
