############################################################################################################################
# Authors: Irene van Blokland, Roy Oelen
# Name: epifat_comparison.R
# Function: to compare with EAT and SAT data of https://academic.oup.com/cardiovascres/article/109/2/228/2196670
############################################################################################################################

###############################
# loading libraries
###############################
library(Seurat)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(dplyr)
library(ggVennDiagram)

####################
# load stuff      
####################
# location of the file
objects_loc <- "/groups/umcg-franke-scrna/tmp01/projects/epifat/ongoing/seurat_preprocess_samples/objects/"
spaceranger_object_loc <- paste(objects_loc, "spaceranger.20230922.integrated.arteries.fulldecon.rds", sep = "")
# read the object
spaceranger_object <- readRDS(spaceranger_object_loc)

####################
# load data      
####################
# external data
compare_eat <- read.table("epifat_compare_EAT.tsv", sep = "\t", header = T)
compare_sat <- read.table("epifat_compare_SAT.tsv", sep = "\t", header = T)

####################
# Venndiagram      
####################

# load EAT and SAT genes from paper
set_A <- compare_endo$Gene.symbol
set_B <- compare_eat$Gene
set_C <- compare_sat$Gene
set_D <-compare_endo$Gene.symbol
# load all our EAT genes
all_st_genes <- rownames(spaceranger_object@assays$SCT@counts)
# generate venndiagram
ggVennDiagram(x = list("eat" = set_B, "sat" = set_C, "all genes" = all_st_genes)) 
