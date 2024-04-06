############################################################################################################################
# Authors: Irene van Blokland
# Name: epifat_semla_nmf.R
# Function: Non-negative matrix factorization according to Semla: https://ludvigla.github.io/semla/articles/NNMF.html
############################################################################################################################

###############################
# loading libraries
###############################
library(ggplot2)
library(dplyr)
library(semla)
library(tibble)
library(singlet)
library(SeuratData)
library(Seurat)

# make sure to use Semla version 1.1.6
install.packages("remotes")
remotes::install_github("ludvigla/semla", ref = "v1.1.6")

###############################
# main code
###############################

# location of the file
objects_loc <- "/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/seurat_preprocess_samples/objects/"
spaceranger_object_loc <- paste(objects_loc, "spaceranger.20230922.integrated.arteries.fulldecon.rds", sep = "")
spaceranger_object_loc <- paste(objects_loc, "spaceranger.20230811.rds", sep = "")

# read the object
spaceranger_object <- readRDS(spaceranger_object_loc)

# make sure we have SCT normalized data
DefaultAssay(spaceranger_object) <- 'SCT'

###############################
# select all sections together
###############################
# for this, the integrated seurat object is taken (spaceranger.20230922.integrated.arteries.fulldecon.rds)
# Update spaceranger_raw for compatibility with semla
# prepare update
spaceranger_semla <- UpdateSeuratForSemla(
  spaceranger_object,
  image_type = c("tissue_lowres"),
  verbose = TRUE
)

# Normalize data and find top variable features
spaceranger_semla <- spaceranger_semla |> 
  NormalizeData() |> 
  FindVariableFeatures()

# OPTIONAL: subset data to improve computational speed
spaceranger_semla <- spaceranger_semla[VariableFeatures(spaceranger_semla), ]

# Set seed for reproducibility
set.seed(42)
# profit
spaceranger_semla <- RunNMF(spaceranger_semla)

# explore the optimal NMF
RankPlot(spaceranger_semla) # 19

# spatial visualization
MapFeatures(spaceranger_semla, 
            features = paste0("NMF_", 1:19), 
            override_plot_dims = TRUE, 
            colors = viridis::magma(n = 11, direction = -1)) &
  theme(plot.title = element_blank())

# map specific features based on eyeballing differences between health and disease
MapFeatures(spaceranger_semla, 
            features = c("NMF_10", "NMF_12", "NMF_17", "NMF_9", "NMF_14", "NMF_19"), 
            override_plot_dims = TRUE, 
            colors = viridis::magma(n = 11, direction = -1)) &
  theme(plot.title = element_blank())

# plot all feature loadings
PlotFeatureLoadings(spaceranger_semla, 
                    dims = 1:19, 
                    reduction = "nmf", 
                    nfeatures = 30,
                    mode = "dotplot", 
                    fill = "darkmagenta",
                    pt_size = 3)

# plot feature loadings for specified NMFs
PlotFeatureLoadings(spaceranger_semla, 
                    dims = c(10, 12, 17, 9, 14, 19), 
                    reduction = "nmf", 
                    nfeatures = 20,
                    mode = "dotplot", 
                    fill = "darkmagenta",
                    pt_size = 3)

spaceranger_semla <- LoadImages(spaceranger_semla)

factor_colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
                   '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                   '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3')

# Select non-overlapping factors
selected_factors <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
# select for NMFs
selected_factors <- c(9, 10, 12, 14, 17, 19)

MapMultipleFeatures(spaceranger_semla, 
                    features = paste0("NMF_", selected_factors), 
                    colors = factor_colors, 
                    image_use = "raw", 
                    override_plot_dims = TRUE, 
                    pt_size = 2)

PlotFeatureLoadings(spaceranger_semla, 
                    dims = selected_factors, 
                    reduction = "nmf", 
                    nfeatures = 10, 
                    mode = "heatmap", 
                    gradient_colors = viridis::magma(n = 11, direction = -1))

# try functional enrichment analysis
# fetch feature.loadings from DimReduc object
nmf_loadings <- spaceranger_semla[["nmf"]]@feature.loadings

# Convert to long format and group data by factor
gene_loadings_sorted <- nmf_loadings |> 
  as.data.frame() |> 
  tibble::rownames_to_column(var = "gene") |> 
  as_tibble() |> 
  tidyr::pivot_longer(all_of(colnames(nmf_loadings)), names_to = "fctr", values_to = "loading") |> 
  mutate(fctr = factor(fctr, colnames(nmf_loadings))) |> 
  group_by(fctr) |> 
  arrange(fctr, -loading)

# Extract top 10 genes per factor
gene_loadings_sorted |> 
  slice_head(n = 10)

library(gprofiler2)

# Get gene sets
gene_set_nmf_1 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_1") |> 
  slice_head(n = 10)
gene_set_nmf_2 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_2") |> 
  slice_head(n = 10)
gene_set_nmf_3 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_3") |> 
  slice_head(n = 10)
gene_set_nmf_4 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_4") |> 
  slice_head(n = 10)
gene_set_nmf_5 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_5") |> 
  slice_head(n = 10)
gene_set_nmf_6 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_6") |> 
  slice_head(n = 10)
gene_set_nmf_7 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_7") |> 
  slice_head(n = 10)
gene_set_nmf_8 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_8") |> 
  slice_head(n = 10)
gene_set_nmf_9 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_9") |> 
  slice_head(n = 10)

gene_set_nmf_10 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_10") |> 
  slice_head(n = 10)
gene_set_nmf_11 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_11") |> 
  slice_head(n = 10)
gene_set_nmf_12 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_12") |> 
  slice_head(n = 10)
gene_set_nmf_13 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_13") |> 
  slice_head(n = 10)
gene_set_nmf_14 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_14") |> 
  slice_head(n = 10)
gene_set_nmf_15 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_15") |> 
  slice_head(n = 10)
gene_set_nmf_16 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_16") |> 
  slice_head(n = 10)
gene_set_nmf_17 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_17") |> 
  slice_head(n = 10)
gene_set_nmf_18 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_18") |> 
  slice_head(n = 10)
gene_set_nmf_19 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_19") |> 
  slice_head(n = 10)

# Get gene sets, get full sets
gene_set_nmf_1 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_1")

gene_set_nmf_2 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_2") 
gene_set_nmf_3 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_3") 
gene_set_nmf_4 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_4")
gene_set_nmf_5 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_5")
gene_set_nmf_6 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_6") 
gene_set_nmf_7 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_7") |> 
  slice_head(n = 10)
gene_set_nmf_8 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_8") |> 
  slice_head(n = 10)
gene_set_nmf_9 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_9") |> 
  slice_head(n = 10)

gene_set_nmf_10 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_10") |> 
  slice_head(n = 10)
gene_set_nmf_11 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_11") |> 
  slice_head(n = 10)
gene_set_nmf_12 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_12") |> 
  slice_head(n = 10)
gene_set_nmf_13 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_13") |> 
  slice_head(n = 10)
gene_set_nmf_14 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_14") |> 
  slice_head(n = 10)
gene_set_nmf_15 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_15") |> 
  slice_head(n = 10)
gene_set_nmf_16 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_16") |> 
  slice_head(n = 10)
gene_set_nmf_17 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_17") |> 
  slice_head(n = 10)
gene_set_nmf_18 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_18") |> 
  slice_head(n = 10)
gene_set_nmf_19 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_19") |> 
  slice_head(n = 10)

# Run FEA
fea_results_nmf_1 <- gost(query = gene_set_nmf_1$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_2 <- gost(query = gene_set_nmf_2$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_3 <- gost(query = gene_set_nmf_3$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_4 <- gost(query = gene_set_nmf_4$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_5 <- gost(query = gene_set_nmf_5$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_6 <- gost(query = gene_set_nmf_6$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_7 <- gost(query = gene_set_nmf_7$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_8 <- gost(query = gene_set_nmf_8$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_9 <- gost(query = gene_set_nmf_9$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_10 <- gost(query = gene_set_nmf_10$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_11 <- gost(query = gene_set_nmf_11$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_12 <- gost(query = gene_set_nmf_12$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_13 <- gost(query = gene_set_nmf_13$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_14 <- gost(query = gene_set_nmf_14$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_15 <- gost(query = gene_set_nmf_15$ene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_16 <- gost(query = gene_set_nmf_16$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_17 <- gost(query = gene_set_nmf_17$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_18 <- gost(query = gene_set_nmf_18$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()
fea_results_nmf_19 <- gost(query = gene_set_nmf_19$gene, ordered_query = TRUE, organism = "hsapiens", sources = "GO:BP")$result |> 
  as_tibble()

# Look at results
fea_results_nmf_1 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_2 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_3 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_4 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_5 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_6 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_7 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_8 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_9 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_10|> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_11 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_12 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_13 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_14 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_15 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_16 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_17 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_18 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_19 |> select(p_value, term_size, query_size, intersection_size, term_name)

###############################
# select section V10A20-016_A1
###############################
# for this, the non-integrated object is taken (spaceranger.20230811.rds)
spaceranger_object_A1 <- spaceranger_object[["V10A20-016_A1"]]

# Update spaceranger_raw for compatibility with semla
# prepare update
spaceranger_semla <- UpdateSeuratForSemla(
  spaceranger_object_A1,
  image_type = c("tissue_lowres"),
  verbose = TRUE
)

# Normalize data and find top variable features
spaceranger_semla <- spaceranger_semla |> 
  NormalizeData() |> 
  FindVariableFeatures()

# OPTIONAL: subset data to improve computational speed
spaceranger_semla <- spaceranger_semla[VariableFeatures(spaceranger_semla), ]

# Set seed for reproducibility
set.seed(42)
# profit
spaceranger_semla <- RunNMF(spaceranger_semla)

# explore the optimal NMF
RankPlot(spaceranger_semla) # 4

# spatial visualization
MapFeatures(spaceranger_semla, 
            features = paste0("NMF_", 1:4), 
            override_plot_dims = TRUE, 
            colors = viridis::magma(n = 11, direction = -1)) &
  theme(plot.title = element_blank())

PlotFeatureLoadings(spaceranger_semla, 
                    dims = 1:2, 
                    reduction = "nmf", 
                    nfeatures = 30,
                    mode = "dotplot", 
                    fill = "darkmagenta",
                    pt_size = 3)

spaceranger_semla <- LoadImages(spaceranger_semla)

factor_colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
                   '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                   '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3')

# Select non-overlapping factors
selected_factors <- c(1, 2, 3, 4)

MapMultipleFeatures(spaceranger_semla, 
                    features = paste0("NMF_", selected_factors), 
                    colors = factor_colors, 
                    image_use = "raw", 
                    override_plot_dims = TRUE, 
                    pt_size = 2)

PlotFeatureLoadings(spaceranger_semla, 
                    dims = selected_factors, 
                    reduction = "nmf", 
                    nfeatures = 30, 
                    mode = "heatmap", 
                    gradient_colors = viridis::magma(n = 11, direction = -1))

# try functional enrichment analysis
# fetch feature.loadings from DimReduc object
nmf_loadings <- spaceranger_semla[["nmf"]]@feature.loadings

# Convert to long format and group data by factor
gene_loadings_sorted <- nmf_loadings |> 
  as.data.frame() |> 
  tibble::rownames_to_column(var = "gene") |> 
  as_tibble() |> 
  tidyr::pivot_longer(all_of(colnames(nmf_loadings)), names_to = "fctr", values_to = "loading") |> 
  mutate(fctr = factor(fctr, colnames(nmf_loadings))) |> 
  group_by(fctr) |> 
  arrange(fctr, -loading)

# Extract top 10 genes per factor
gene_loadings_sorted |> 
  slice_head(n = 10)

library(gprofiler2)

# Get gene sets
gene_set_nmf_1 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_1") |> 
  slice_head(n = 10)
gene_set_nmf_2 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_2") |> 
  slice_head(n = 10)
gene_set_nmf_3 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_3") |> 
  slice_head(n = 10)
gene_set_nmf_4 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_4") |> 
  slice_head(n = 10)

# Run FEA
fea_results_nmf_1 <- gost(query = gene_set_nmf_1$gene, ordered_query = TRUE, organism = "hsapiens", sources = "REAC")$result |> 
  as_tibble()
fea_results_nmf_2 <- gost(query = gene_set_nmf_2$gene, ordered_query = TRUE, organism = "hsapiens", sources = "REAC")$result |> 
  as_tibble()
fea_results_nmf_3 <- gost(query = gene_set_nmf_3$gene, ordered_query = TRUE, organism = "hsapiens", sources = "REAC")$result |> 
  as_tibble()
fea_results_nmf_4 <- gost(query = gene_set_nmf_4$gene, ordered_query = TRUE, organism = "hsapiens", sources = "REAC")$result |> 
  as_tibble()

# Look at results
fea_results_nmf_1 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_2 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_3 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_4 |> select(p_value, term_size, query_size, intersection_size, term_name)


###############################
# select section V10A20-016_B1
###############################
spaceranger_object_B1 <- spaceranger_object[["V10A20-016_B1"]]

# Update spaceranger_raw for compatibility with semla
# prepare update
spaceranger_semla <- UpdateSeuratForSemla(
  spaceranger_object_B1,
  image_type = c("tissue_lowres"),
  verbose = TRUE
)

# Normalize data and find top variable features
spaceranger_semla <- spaceranger_semla |> 
  NormalizeData() |> 
  FindVariableFeatures()

# OPTIONAL: subset data to improve computational speed
spaceranger_semla <- spaceranger_semla[VariableFeatures(spaceranger_semla), ]

# Set seed for reproducibility
set.seed(42)
# profit
spaceranger_semla <- RunNMF(spaceranger_semla)

# explore the optimal NMF
RankPlot(spaceranger_semla) # 6

# spatial visualization
MapFeatures(spaceranger_semla, 
            features = paste0("NMF_", 1:6), 
            override_plot_dims = TRUE, 
            colors = viridis::magma(n = 11, direction = -1)) &
  theme(plot.title = element_blank())

PlotFeatureLoadings(spaceranger_semla, 
                    dims = 1:2, 
                    reduction = "nmf", 
                    nfeatures = 30,
                    mode = "dotplot", 
                    fill = "darkmagenta",
                    pt_size = 3)

spaceranger_semla <- LoadImages(spaceranger_semla)

factor_colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
                   '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                   '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3')

# Select non-overlapping factors
selected_factors <- c(1, 2, 3, 4, 5, 6)

MapMultipleFeatures(spaceranger_semla, 
                    features = paste0("NMF_", selected_factors), 
                    colors = factor_colors, 
                    image_use = "raw", 
                    override_plot_dims = TRUE, 
                    pt_size = 2)

PlotFeatureLoadings(spaceranger_semla, 
                    dims = selected_factors, 
                    reduction = "nmf", 
                    nfeatures = 20, 
                    mode = "heatmap", 
                    gradient_colors = viridis::magma(n = 11, direction = -1))

# try functional enrichment analysis
# fetch feature.loadings from DimReduc object
nmf_loadings <- spaceranger_semla[["nmf"]]@feature.loadings

# Convert to long format and group data by factor
gene_loadings_sorted <- nmf_loadings |> 
  as.data.frame() |> 
  tibble::rownames_to_column(var = "gene") |> 
  as_tibble() |> 
  tidyr::pivot_longer(all_of(colnames(nmf_loadings)), names_to = "fctr", values_to = "loading") |> 
  mutate(fctr = factor(fctr, colnames(nmf_loadings))) |> 
  group_by(fctr) |> 
  arrange(fctr, -loading)

# Extract top 10 genes per factor
gene_loadings_sorted |> 
  slice_head(n = 10)

library(gprofiler2)

# Get gene sets
gene_set_nmf_1 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_1") |> 
  slice_head(n = 10)
gene_set_nmf_2 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_2") |> 
  slice_head(n = 10)
gene_set_nmf_3 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_3") |> 
  slice_head(n = 10)
gene_set_nmf_4 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_4") |> 
  slice_head(n = 10)
gene_set_nmf_5 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_5") |> 
  slice_head(n = 10)
gene_set_nmf_6 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_6") |> 
  slice_head(n = 10)

# Run FEA
fea_results_nmf_1 <- gost(query = gene_set_nmf_1$gene, ordered_query = TRUE, organism = "hsapiens", sources = "REAC")$result |> 
  as_tibble()
fea_results_nmf_2 <- gost(query = gene_set_nmf_2$gene, ordered_query = TRUE, organism = "hsapiens", sources = "REAC")$result |> 
  as_tibble()
fea_results_nmf_3 <- gost(query = gene_set_nmf_3$gene, ordered_query = TRUE, organism = "hsapiens", sources = "REAC")$result |> 
  as_tibble()
fea_results_nmf_4 <- gost(query = gene_set_nmf_4$gene, ordered_query = TRUE, organism = "hsapiens", sources = "REAC")$result |> 
  as_tibble()
fea_results_nmf_5 <- gost(query = gene_set_nmf_5$gene, ordered_query = TRUE, organism = "hsapiens", sources = "REAC")$result |> 
  as_tibble()
fea_results_nmf_6 <- gost(query = gene_set_nmf_6$gene, ordered_query = TRUE, organism = "hsapiens", sources = "REAC")$result |> 
  as_tibble()


# Look at results
fea_results_nmf_1 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_2 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_3 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_4 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_5 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_6 |> select(p_value, term_size, query_size, intersection_size, term_name)


###############################
# select section V10A20-016_C1
###############################
spaceranger_object_C1 <- spaceranger_object[["V10A20-016_C1"]]

# Update spaceranger_raw for compatibility with semla
# prepare update
spaceranger_semla <- UpdateSeuratForSemla(
  spaceranger_object_C1,
  image_type = c("tissue_lowres"),
  verbose = TRUE
)

# Normalize data and find top variable features
spaceranger_semla <- spaceranger_semla |> 
  NormalizeData() |> 
  FindVariableFeatures()

# OPTIONAL: subset data to improve computational speed
spaceranger_semla <- spaceranger_semla[VariableFeatures(spaceranger_semla), ]

# Set seed for reproducibility
set.seed(42)
# profit
spaceranger_semla <- RunNMF(spaceranger_semla)

# explore the optimal NMF
RankPlot(spaceranger_semla) # 7

# spatial visualization
MapFeatures(spaceranger_semla, 
            features = paste0("NMF_", 1:7), 
            override_plot_dims = TRUE, 
            colors = viridis::magma(n = 11, direction = -1)) &
  theme(plot.title = element_blank())

PlotFeatureLoadings(spaceranger_semla, 
                    dims = 1:2, 
                    reduction = "nmf", 
                    nfeatures = 30,
                    mode = "dotplot", 
                    fill = "darkmagenta",
                    pt_size = 3)

spaceranger_semla <- LoadImages(spaceranger_semla)

factor_colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
                   '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                   '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3')

# Select non-overlapping factors
selected_factors <- c(1, 2, 3, 4, 5, 6, 7)

MapMultipleFeatures(spaceranger_semla, 
                    features = paste0("NMF_", selected_factors), 
                    colors = factor_colors, 
                    image_use = "raw", 
                    override_plot_dims = TRUE, 
                    pt_size = 2)

PlotFeatureLoadings(spaceranger_semla, 
                    dims = selected_factors, 
                    reduction = "nmf", 
                    nfeatures = 10, 
                    mode = "heatmap", 
                    gradient_colors = viridis::magma(n = 11, direction = -1))

# try functional enrichment analysis
# fetch feature.loadings from DimReduc object
nmf_loadings <- spaceranger_semla[["nmf"]]@feature.loadings

# Convert to long format and group data by factor
gene_loadings_sorted <- nmf_loadings |> 
  as.data.frame() |> 
  tibble::rownames_to_column(var = "gene") |> 
  as_tibble() |> 
  tidyr::pivot_longer(all_of(colnames(nmf_loadings)), names_to = "fctr", values_to = "loading") |> 
  mutate(fctr = factor(fctr, colnames(nmf_loadings))) |> 
  group_by(fctr) |> 
  arrange(fctr, -loading)

# Extract top 10 genes per factor
gene_loadings_sorted |> 
  slice_head(n = 10)

library(gprofiler2)

# Get gene sets
gene_set_nmf_1 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_1") |> 
  slice_head(n = 10)
gene_set_nmf_2 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_2") |> 
  slice_head(n = 10)
gene_set_nmf_3 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_3") |> 
  slice_head(n = 10)
gene_set_nmf_4 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_4") |> 
  slice_head(n = 10)
gene_set_nmf_5 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_5") |> 
  slice_head(n = 10)
gene_set_nmf_6 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_6") |> 
  slice_head(n = 10)
gene_set_nmf_7 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_7") |> 
  slice_head(n = 10)
# Run FEA
fea_results_nmf_1 <- gost(query = gene_set_nmf_1$gene, ordered_query = TRUE, organism = "hsapiens", sources = "KEGG")$result |> 
  as_tibble()
fea_results_nmf_2 <- gost(query = gene_set_nmf_2$gene, ordered_query = TRUE, organism = "hsapiens", sources = "KEGG")$result |> 
  as_tibble()
fea_results_nmf_3 <- gost(query = gene_set_nmf_3$gene, ordered_query = TRUE, organism = "hsapiens", sources = "KEGG")$result |> 
  as_tibble()
fea_results_nmf_4 <- gost(query = gene_set_nmf_4$gene, ordered_query = TRUE, organism = "hsapiens", sources = "KEGG")$result |> 
  as_tibble()
fea_results_nmf_5 <- gost(query = gene_set_nmf_5$gene, ordered_query = TRUE, organism = "hsapiens", sources = "KEGG")$result |> 
  as_tibble()
fea_results_nmf_6 <- gost(query = gene_set_nmf_6$gene, ordered_query = TRUE, organism = "hsapiens", sources = "KEGG")$result |> 
  as_tibble()
fea_results_nmf_7 <- gost(query = gene_set_nmf_7$gene, ordered_query = TRUE, organism = "hsapiens", sources = "KEGG")$result |> 
  as_tibble()

# Look at results
fea_results_nmf_1 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_2 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_3 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_4 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_5 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_6 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_7 |> select(p_value, term_size, query_size, intersection_size, term_name)

###############################
# select section V10A20-016_D1
###############################
spaceranger_object_D1 <- spaceranger_object[["V10A20-016_D1"]]

# Update spaceranger_raw for compatibility with semla
# prepare update
spaceranger_semla <- UpdateSeuratForSemla(
  spaceranger_object_D1,
  image_type = c("tissue_lowres"),
  verbose = TRUE
)

# Normalize data and find top variable features
spaceranger_semla <- spaceranger_semla |> 
  NormalizeData() |> 
  FindVariableFeatures()

# OPTIONAL: subset data to improve computational speed
spaceranger_semla <- spaceranger_semla[VariableFeatures(spaceranger_semla), ]

# Set seed for reproducibility
set.seed(42)
# profit
spaceranger_semla <- RunNMF(spaceranger_semla)

# explore the optimal NMF
RankPlot(spaceranger_semla) # 6

# spatial visualization
MapFeatures(spaceranger_semla, 
            features = paste0("NMF_", 1:6), 
            override_plot_dims = TRUE, 
            colors = viridis::magma(n = 11, direction = -1)) &
  theme(plot.title = element_blank())

PlotFeatureLoadings(spaceranger_semla, 
                    dims = 1:2, 
                    reduction = "nmf", 
                    nfeatures = 30,
                    mode = "dotplot", 
                    fill = "darkmagenta",
                    pt_size = 3)

spaceranger_semla <- LoadImages(spaceranger_semla)

factor_colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
                   '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                   '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3')

# Select non-overlapping factors
selected_factors <- c(1, 2, 3, 4, 5, 6)

MapMultipleFeatures(spaceranger_semla, 
                    features = paste0("NMF_", selected_factors), 
                    colors = factor_colors, 
                    image_use = "raw", 
                    override_plot_dims = TRUE, 
                    pt_size = 2)

PlotFeatureLoadings(spaceranger_semla, 
                    dims = selected_factors, 
                    reduction = "nmf", 
                    nfeatures = 20, 
                    mode = "heatmap", 
                    gradient_colors = viridis::magma(n = 11, direction = -1))

# try functional enrichment analysis
# fetch feature.loadings from DimReduc object
nmf_loadings <- spaceranger_semla[["nmf"]]@feature.loadings

# Convert to long format and group data by factor
gene_loadings_sorted <- nmf_loadings |> 
  as.data.frame() |> 
  tibble::rownames_to_column(var = "gene") |> 
  as_tibble() |> 
  tidyr::pivot_longer(all_of(colnames(nmf_loadings)), names_to = "fctr", values_to = "loading") |> 
  mutate(fctr = factor(fctr, colnames(nmf_loadings))) |> 
  group_by(fctr) |> 
  arrange(fctr, -loading)

# Extract top 10 genes per factor
gene_loadings_sorted |> 
  slice_head(n = 10)

library(gprofiler2)

# Get gene sets
gene_set_nmf_1 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_1") |> 
  slice_head(n = 10)
gene_set_nmf_2 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_2") |> 
  slice_head(n = 10)
gene_set_nmf_3 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_3") |> 
  slice_head(n = 10)
gene_set_nmf_4 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_4") |> 
  slice_head(n = 10)
gene_set_nmf_5 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_5") |> 
  slice_head(n = 10)
gene_set_nmf_6 <- gene_loadings_sorted |> 
  filter(fctr == "NMF_6") |> 
  slice_head(n = 10)

# Run FEA
fea_results_nmf_1 <- gost(query = gene_set_nmf_1$gene, ordered_query = TRUE, organism = "hsapiens", sources = "KEGG")$result |> 
  as_tibble()
fea_results_nmf_2 <- gost(query = gene_set_nmf_2$gene, ordered_query = TRUE, organism = "hsapiens", sources = "KEGG")$result |> 
  as_tibble()
fea_results_nmf_3 <- gost(query = gene_set_nmf_3$gene, ordered_query = TRUE, organism = "hsapiens", sources = "KEGG")$result |> 
  as_tibble()
fea_results_nmf_4 <- gost(query = gene_set_nmf_4$gene, ordered_query = TRUE, organism = "hsapiens", sources = "KEGG")$result |> 
  as_tibble()
fea_results_nmf_5 <- gost(query = gene_set_nmf_5$gene, ordered_query = TRUE, organism = "hsapiens", sources = "KEGG")$result |> 
  as_tibble()
fea_results_nmf_6 <- gost(query = gene_set_nmf_6$gene, ordered_query = TRUE, organism = "hsapiens", sources = "KEGG")$result |> 
  as_tibble()

# Look at results
fea_results_nmf_1 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_2 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_3 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_4 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_5 |> select(p_value, term_size, query_size, intersection_size, term_name)
fea_results_nmf_6 |> select(p_value, term_size, query_size, intersection_size, term_name)

