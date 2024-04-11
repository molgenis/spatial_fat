# script by Oelen and van Blokland

###############################
# load libraries
###############################

library(Seurat) #"should be Seurat version 5!"
library(Matrix)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggpubr)

###############################
# functions
###############################

get_decon_data_per_section_and_cluster <- function(metadata, cell_types, celltype_column_prepend="", section_column="section", cluster_column="seurat_clusters") {
  # store per section
  data_per_section <- list()
  # check each section
  for (section in unique(metadata[[section_column]])) {
    # subset to section
    metadata_section <- metadata[!is.na(metadata[[section_column]]) & metadata[[section_column]] == section, ]
    # we will store for every cluster
    data_per_cluster <- list()
    # check each cluster
    for (cluster in unique(metadata_section[[cluster_column]])) {
      # subset to cluster
      metadata_cluster <- metadata_section[!is.na(metadata_section[[cluster_column]]) & metadata_section[[cluster_column]] == cluster, ]
      # check each cell type and put results in list
      celltype_list <- list()
      expression_list <- list()
      for (cell_type in cell_types) {
        # add a dataframe with the cell type names
        celltype_list[[cell_type]] <- data.frame(cell_type = rep(cell_type, times = nrow(metadata_cluster)))
        # make the actual column name
        celltype_column_name <- paste(celltype_column_prepend, cell_type, sep = "")
        # get the value for that column, and put into dataframe
        expression_list[[cell_type]] <- data.frame(expression = metadata_cluster[[celltype_column_name]])
      }
      # merge all cell types together
      celltype_expressions <- cbind(
        do.call("rbind", celltype_list),
        do.call("rbind", expression_list)
      )
      # put the result in the data for the cluster
      data_per_cluster[[as.character(cluster)]] <- celltype_expressions
    }
    # put the results in the data for the section
    data_per_section[[section]] <- data_per_cluster
  }
  return(data_per_section)
}


plot_expression_density <- function(input_table, cell_type_column="cell_type", expression_column="expression") {
  # Create density plot using ggplot2
  p <- ggplot(data = NULL, aes(x = input_table[[expression_column]], fill = input_table[[cell_type_column]])) +
    geom_density(alpha = 0.7) +
    ggtitle("Density Plot") +
    xlab("Expression") +
    ylab("Density")
  return(p)
}

###############################
# main
###############################

# read the object 
spaceranger_tangram <- readRDS("/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/seurat_preprocess_samples/objects/spaceranger.20230922.integrated.arteries.fulldecon.rds")

#########################
# plotting density plots
#########################
# Select cell types
selected_cell_types <- c("Adipocyte", "Atrial_Cardiomyocyte", "Endothelial_cell", 
                         "Fibroblast", "Lymphatic_Endothelial_cell", "Lymphoid", 
                         "Mast_cell", "Mesothelial_cell", "Mural_cell", 
                         "Myeloid", "Neural_cell", "Ventricular_Cardiomyocyte")

# get the decon data per cluster and section for RCTD and Tangram
rctd_per_section_cluster <- get_decon_data_per_section_and_cluster(spaceranger_tangram@meta.data, selected_cell_types)
tangram_per_section_cluster <- get_decon_data_per_section_and_cluster(spaceranger_tangram@meta.data, selected_cell_types, celltype_column_prepend = "tangram.")

# plot for RCTD per section
plot_grid(
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_A1"]][["0"]]) + ggtitle("V10A20_016_A1, cluster 0") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_A1"]][["1"]]) + ggtitle("V10A20_016_A1, cluster 1") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_A1"]][["2"]]) + ggtitle("V10A20_016_A1, cluster 2") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_A1"]][["3"]]) + ggtitle("V10A20_016_A1, cluster 3") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_A1"]][["4"]]) + ggtitle("V10A20_016_A1, cluster 4") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_A1"]][["5"]]) + ggtitle("V10A20_016_A1, cluster 5") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_A1"]][["6"]]) + ggtitle("V10A20_016_A1, cluster 6") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_A1"]][["7"]]) + ggtitle("V10A20_016_A1, cluster 7") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_A1"]][["8"]]) + ggtitle("V10A20_016_A1, cluster 8") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_A1"]][["9"]]) + ggtitle("V10A20_016_A1, cluster 9") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_A1"]][["10"]]) + ggtitle("V10A20_016_A1, cluster 10") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_A1"]][["11"]]) + ggtitle("V10A20_016_A1, cluster 11") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_A1"]][["12"]]) + ggtitle("V10A20_016_A1, cluster 12") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_A1"]][["13"]]) + ggtitle("V10A20_016_A1, cluster 13") + theme(legend.position = "none") + ylim(0, 25),
  get_legend(plot_expression_density(rctd_per_section_cluster[["V10A20_016_A1"]][["13"]]))
)

plot_grid(
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_B1"]][["0"]]) + ggtitle("V10A20_016_B1, cluster 0") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_B1"]][["1"]]) + ggtitle("V10A20_016_B1, cluster 1") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_B1"]][["2"]]) + ggtitle("V10A20_016_B1, cluster 2") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_B1"]][["3"]]) + ggtitle("V10A20_016_B1, cluster 3") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_B1"]][["4"]]) + ggtitle("V10A20_016_B1, cluster 4") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_B1"]][["5"]]) + ggtitle("V10A20_016_B1, cluster 5") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_B1"]][["6"]]) + ggtitle("V10A20_016_B1, cluster 6") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_B1"]][["7"]]) + ggtitle("V10A20_016_B1, cluster 7") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_B1"]][["8"]]) + ggtitle("V10A20_016_B1, cluster 8") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_B1"]][["9"]]) + ggtitle("V10A20_016_B1, cluster 9") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_B1"]][["10"]]) + ggtitle("V10A20_016_B1, cluster 10") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_B1"]][["11"]]) + ggtitle("V10A20_016_B1, cluster 11") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_B1"]][["12"]]) + ggtitle("V10A20_016_B1, cluster 12") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_B1"]][["13"]]) + ggtitle("V10A20_016_B1, cluster 13") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_B1"]][["14"]]) + ggtitle("V10A20_016_B1, cluster 14") + theme(legend.position = "none") + ylim(0, 25),
  get_legend(plot_expression_density(rctd_per_section_cluster[["V10A20_016_B1"]][["14"]]))
)

plot_grid(
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_C1"]][["0"]]) + ggtitle("V10A20_016_C1, cluster 0") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_C1"]][["1"]]) + ggtitle("V10A20_016_C1, cluster 1") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_C1"]][["2"]]) + ggtitle("V10A20_016_C1, cluster 2") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_C1"]][["3"]]) + ggtitle("V10A20_016_C1, cluster 3") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_C1"]][["4"]]) + ggtitle("V10A20_016_C1, cluster 4") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_C1"]][["5"]]) + ggtitle("V10A20_016_C1, cluster 5") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_C1"]][["6"]]) + ggtitle("V10A20_016_C1, cluster 6") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_C1"]][["7"]]) + ggtitle("V10A20_016_C1, cluster 7") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_C1"]][["8"]]) + ggtitle("V10A20_016_C1, cluster 8") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_C1"]][["9"]]) + ggtitle("V10A20_016_C1, cluster 9") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_C1"]][["10"]]) + ggtitle("V10A20_016_C1, cluster 10") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_C1"]][["11"]]) + ggtitle("V10A20_016_C1, cluster 11") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_C1"]][["12"]]) + ggtitle("V10A20_016_C1, cluster 12") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_C1"]][["13"]]) + ggtitle("V10A20_016_C1, cluster 13") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_C1"]][["14"]]) + ggtitle("V10A20_016_C1, cluster 14") + theme(legend.position = "none") + ylim(0, 25),
  get_legend(plot_expression_density(rctd_per_section_cluster[["V10A20_016_C1"]][["14"]]))
)

plot_grid(
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_D1"]][["0"]]) + ggtitle("V10A20_016_D1, cluster 0") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_D1"]][["1"]]) + ggtitle("V10A20_016_D1, cluster 1") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_D1"]][["2"]]) + ggtitle("V10A20_016_D1, cluster 2") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_D1"]][["3"]]) + ggtitle("V10A20_016_D1, cluster 3") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_D1"]][["4"]]) + ggtitle("V10A20_016_D1, cluster 4") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_D1"]][["5"]]) + ggtitle("V10A20_016_D1, cluster 5") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_D1"]][["6"]]) + ggtitle("V10A20_016_D1, cluster 6") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_D1"]][["7"]]) + ggtitle("V10A20_016_D1, cluster 7") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_D1"]][["8"]]) + ggtitle("V10A20_016_D1, cluster 8") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_D1"]][["9"]]) + ggtitle("V10A20_016_D1, cluster 9") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_D1"]][["10"]]) + ggtitle("V10A20_016_D1, cluster 10") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_D1"]][["11"]]) + ggtitle("V10A20_016_D1, cluster 11") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_D1"]][["12"]]) + ggtitle("V10A20_016_D1, cluster 12") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_D1"]][["13"]]) + ggtitle("V10A20_016_D1, cluster 13") + theme(legend.position = "none") + ylim(0, 25),
  plot_expression_density(rctd_per_section_cluster[["V10A20_016_D1"]][["14"]]) + ggtitle("V10A20_016_D1, cluster 14") + theme(legend.position = "none") + ylim(0, 25),
  get_legend(plot_expression_density(rctd_per_section_cluster[["V10A20_016_D1"]][["14"]]))
)

# plot for tangram per section
plot_grid(
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_A1"]][["0"]]) + ggtitle("V10A20_016_A1, cluster 0") + theme(legend.position = "none") + ylim(0, 0.1) + xlim(0, 125),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_A1"]][["1"]]) + ggtitle("V10A20_016_A1, cluster 1") + theme(legend.position = "none") + ylim(0, .75) + xlim(0, 40),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_A1"]][["2"]]) + ggtitle("V10A20_016_A1, cluster 2") + theme(legend.position = "none") + ylim(0, .4) + xlim(0, 80),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_A1"]][["3"]]) + ggtitle("V10A20_016_A1, cluster 3") + theme(legend.position = "none") + ylim(0, .2) + xlim(0, 75),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_A1"]][["4"]]) + ggtitle("V10A20_016_A1, cluster 4") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_A1"]][["5"]]) + ggtitle("V10A20_016_A1, cluster 5") + theme(legend.position = "none") + ylim(0, .3) + xlim(0, 80),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_A1"]][["6"]]) + ggtitle("V10A20_016_A1, cluster 6") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_A1"]][["7"]]) + ggtitle("V10A20_016_A1, cluster 7") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 100),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_A1"]][["8"]]) + ggtitle("V10A20_016_A1, cluster 8") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 100),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_A1"]][["9"]]) + ggtitle("V10A20_016_A1, cluster 9") + theme(legend.position = "none") + ylim(0, .4) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_A1"]][["10"]]) + ggtitle("V10A20_016_A1, cluster 10") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_A1"]][["11"]]) + ggtitle("V10A20_016_A1, cluster 11") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_A1"]][["12"]]) + ggtitle("V10A20_016_A1, cluster 12") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 80),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_A1"]][["13"]]) + ggtitle("V10A20_016_A1, cluster 13") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 80),
  get_legend(plot_expression_density(tangram_per_section_cluster[["V10A20_016_A1"]][["13"]]))
)

plot_grid(
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_B1"]][["0"]]) + ggtitle("V10A20_016_B1, cluster 0") + theme(legend.position = "none") + ylim(0, 0.1) + xlim(0, 125),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_B1"]][["1"]]) + ggtitle("V10A20_016_B1, cluster 1") + theme(legend.position = "none") + ylim(0, .75) + xlim(0, 40),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_B1"]][["2"]]) + ggtitle("V10A20_016_B1, cluster 2") + theme(legend.position = "none") + ylim(0, .4) + xlim(0, 80),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_B1"]][["3"]]) + ggtitle("V10A20_016_B1, cluster 3") + theme(legend.position = "none") + ylim(0, .2) + xlim(0, 75),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_B1"]][["4"]]) + ggtitle("V10A20_016_B1, cluster 4") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_B1"]][["5"]]) + ggtitle("V10A20_016_B1, cluster 5") + theme(legend.position = "none") + ylim(0, .3) + xlim(0, 80),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_B1"]][["6"]]) + ggtitle("V10A20_016_B1, cluster 6") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_B1"]][["7"]]) + ggtitle("V10A20_016_B1, cluster 7") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 100),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_B1"]][["8"]]) + ggtitle("V10A20_016_B1, cluster 8") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 100),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_B1"]][["9"]]) + ggtitle("V10A20_016_B1, cluster 9") + theme(legend.position = "none") + ylim(0, .4) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_B1"]][["10"]]) + ggtitle("V10A20_016_B1, cluster 10") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_B1"]][["11"]]) + ggtitle("V10A20_016_B1, cluster 11") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_B1"]][["12"]]) + ggtitle("V10A20_016_B1, cluster 12") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 80),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_B1"]][["13"]]) + ggtitle("V10A20_016_B1, cluster 13") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 80),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_B1"]][["14"]]) + ggtitle("V10A20_016_B1, cluster 14") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 80),
  get_legend(plot_expression_density(tangram_per_section_cluster[["V10A20_016_B1"]][["13"]]))
)

plot_grid(
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_C1"]][["0"]]) + ggtitle("V10A20_016_C1, cluster 0") + theme(legend.position = "none") + ylim(0, 0.1) + xlim(0, 125),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_C1"]][["1"]]) + ggtitle("V10A20_016_C1, cluster 1") + theme(legend.position = "none") + ylim(0, .75) + xlim(0, 40),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_C1"]][["2"]]) + ggtitle("V10A20_016_C1, cluster 2") + theme(legend.position = "none") + ylim(0, .4) + xlim(0, 80),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_C1"]][["3"]]) + ggtitle("V10A20_016_C1, cluster 3") + theme(legend.position = "none") + ylim(0, .2) + xlim(0, 75),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_C1"]][["4"]]) + ggtitle("V10A20_016_C1, cluster 4") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_C1"]][["5"]]) + ggtitle("V10A20_016_C1, cluster 5") + theme(legend.position = "none") + ylim(0, .3) + xlim(0, 80),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_C1"]][["6"]]) + ggtitle("V10A20_016_C1, cluster 6") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_C1"]][["7"]]) + ggtitle("V10A20_016_C1, cluster 7") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 100),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_C1"]][["8"]]) + ggtitle("V10A20_016_C1, cluster 8") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 100),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_C1"]][["9"]]) + ggtitle("V10A20_016_C1, cluster 9") + theme(legend.position = "none") + ylim(0, .4) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_C1"]][["10"]]) + ggtitle("V10A20_016_C1, cluster 10") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_C1"]][["11"]]) + ggtitle("V10A20_016_C1, cluster 11") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_C1"]][["12"]]) + ggtitle("V10A20_016_C1, cluster 12") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 80),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_C1"]][["13"]]) + ggtitle("V10A20_016_C1, cluster 13") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 80),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_C1"]][["14"]]) + ggtitle("V10A20_016_C1, cluster 14") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 80),
  get_legend(plot_expression_density(tangram_per_section_cluster[["V10A20_016_C1"]][["13"]]))
)

plot_grid(
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_D1"]][["0"]]) + ggtitle("V10A20_016_D1, cluster 0") + theme(legend.position = "none") + ylim(0, 0.1) + xlim(0, 125),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_D1"]][["1"]]) + ggtitle("V10A20_016_D1, cluster 1") + theme(legend.position = "none") + ylim(0, .75) + xlim(0, 40),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_D1"]][["2"]]) + ggtitle("V10A20_016_D1, cluster 2") + theme(legend.position = "none") + ylim(0, .4) + xlim(0, 80),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_D1"]][["3"]]) + ggtitle("V10A20_016_D1, cluster 3") + theme(legend.position = "none") + ylim(0, .2) + xlim(0, 75),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_D1"]][["4"]]) + ggtitle("V10A20_016_D1, cluster 4") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_D1"]][["5"]]) + ggtitle("V10A20_016_D1, cluster 5") + theme(legend.position = "none") + ylim(0, .3) + xlim(0, 80),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_D1"]][["6"]]) + ggtitle("V10A20_016_D1, cluster 6") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_D1"]][["7"]]) + ggtitle("V10A20_016_D1, cluster 7") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 100),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_D1"]][["8"]]) + ggtitle("V10A20_016_D1, cluster 8") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 100),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_D1"]][["9"]]) + ggtitle("V10A20_016_D1, cluster 9") + theme(legend.position = "none") + ylim(0, .4) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_D1"]][["10"]]) + ggtitle("V10A20_016_D1, cluster 10") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_D1"]][["11"]]) + ggtitle("V10A20_016_D1, cluster 11") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 50),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_D1"]][["12"]]) + ggtitle("V10A20_016_D1, cluster 12") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 80),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_D1"]][["13"]]) + ggtitle("V10A20_016_D1, cluster 13") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 80),
  plot_expression_density(tangram_per_section_cluster[["V10A20_016_D1"]][["14"]]) + ggtitle("V10A20_016_D1, cluster 14") + theme(legend.position = "none") + ylim(0, .25) + xlim(0, 80),
  get_legend(plot_expression_density(tangram_per_section_cluster[["V10A20_016_D1"]][["13"]]))
)




