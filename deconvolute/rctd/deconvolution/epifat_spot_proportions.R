############################################################################################################################
# Authors: Irene van Blokland, Roy Oelen
# Name: epifat_spotprop_cluster.R
# Function: To make a plot x=section/state, y= spot proportion assigned to clusters (%)
############################################################################################################################

###############################
# loading libraries
###############################
library(Seurat)
library(Matrix)
library(ggplot2)

###############################
# loading objects
###############################
# location of the file
objects_loc <- “/groups/umcg-franke-scrna/tmp02/projects/epifat/ongoing/seurat_preprocess_samples/objects/”
spaceranger_integrated_loc <- paste(objects_loc, “spaceranger.20230816.integrated.arteries.rds”, sep = “”)
# read the object
spaceranger_integrated <- readRDS(spaceranger_integrated_loc)

###############################
# loading generating plot
###############################
# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
# Get the cluster assignment from Seurat clusters
seurat_clusters <- spaceranger_integrated@meta.data$seurat_clusters
# Create a data frame combining “section,” “status,” and “seurat_clusters”
data_df <- spaceranger_integrated@meta.data %>%
  select(section, status, seurat_clusters)
# Create a data frame combining “sections A1, C1,” “status,” and “seurat_clusters”
data_df_status <- spaceranger_integrated@meta.data %>%
  filter(section %in% c(“V10A20_016_A1", “V10A20_016_C1”)) %>%
  select(section, status, seurat_clusters)
# Create a data frame combining all sections and “status,” and “seurat_clusters”
data_df_status_all <- spaceranger_integrated@meta.data %>%
  filter(section %in% c(“V10A20_016_A1", “V10A20_016_B1”, “V10A20_016_C1", “V10A20_016_D1”)) %>%
  select(section, status, seurat_clusters)
# Get all unique clusters
all_clusters <- unique(seurat_clusters)
# Calculate the number of spots (cells) assigned to each cluster for each combination of section and status
cluster_counts <- data_df %>%
  group_by(section, seurat_clusters) %>%
  summarize(cluster_count = sum(seurat_clusters %in% all_clusters)) %>%
  ungroup()  # Remove grouping
cluster_counts <- data_df_status %>%
  group_by(section, status, seurat_clusters) %>%
  summarize(cluster_count = sum(seurat_clusters %in% all_clusters)) %>%
  ungroup()  # Remove grouping
cluster_counts <- data_df_status_all %>%
  group_by(section, status, seurat_clusters) %>%
  summarize(cluster_count = sum(seurat_clusters %in% all_clusters)) %>%
  ungroup()  # Remove grouping
# Create a stacked barplot with spot counts using ggplot2 and the custom palette
ggplot(cluster_counts, aes(x = interaction(section, status), y = cluster_count, fill = factor(seurat_clusters))) +
  geom_bar(stat = “identity”, position = position_fill(reverse = TRUE)) +
  labs(x = “Section & Status”, y = “Cluster Count”) +
  scale_fill_manual(values = c(“#EE766F”, “#E58702", “#C99904”, “#A3A51C”, “#6CB22D”, “#3DAB3D”, “#3BB07A”, “#34B5A9", “#1BB8D5”,“#38ABE2",“#6E94CD”, “#9984BD”, “#B377B1”, “#D26EA8", “#EC69A3”)) +  # Use the custom color palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(“Cluster Counts by Section and Status”)
ggplot(cluster_counts, aes(x = interaction(section), y = cluster_count, fill = factor(seurat_clusters))) +
  geom_bar(stat = “identity”, position = position_fill(reverse = TRUE)) +
  labs(x = “Section & Status”, y = “Cluster Count”) +
  scale_fill_manual(values = c(“#EE766F”, “#E58702”, “#C99904", “#A3A51C”, “#6CB22D”, “#3DAB3D”, “#3BB07A”, “#34B5A9”, “#1BB8D5",“#38ABE2”,“#6E94CD”, “#9984BD”, “#B377B1", “#D26EA8”, “#EC69A3")) +  # Use the custom color palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(“Cluster Counts by Section and Status”)
# Calculate the number of spots (cells) assigned to each cluster for each combination of section and status
cluster_counts <- data_df_status %>%
  group_by(section, status, seurat_clusters) %>%
  summarize(cluster_count = sum(seurat_clusters %in% all_clusters)) %>%
  ungroup()  # Remove grouping
# Create a stacked barplot with spot counts using ggplot2 and the custom palette
ggplot(cluster_counts, aes(x = interaction(section, status, drop = TRUE), y = cluster_count, fill = factor(seurat_clusters))) +
  geom_bar(stat = “identity”, position = “stack”) +
  labs(x = “Section & Status”, y = “Cluster Count”) +
  scale_fill_manual(values = c(“#EE766F”, “#E58702”, “#C99904", “#A3A51C”, “#6CB22D”, “#3DAB3D”, “#3BB07A”, “#34B5A9”, “#1BB8D5", “#38ABE2”, “#6E94CD”, “#9984BD”, “#B377B1", “#D26EA8”, “#EC69A3")) +  # Use the custom color palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(“Cluster Counts by Section and Status”)
  
###########################
# Figure 2B: Stacked plot of cell types
###########################
library(ggplot2)
library(tidyr)
library(dplyr)
# Create a stacked barplot with percentage cells per section 
seurat_df <- spaceranger_integrated@meta.data %>%
  select(section, status, Adipocyte, Atrial_Cardiomyocyte, Endothelial_cell, Fibroblast,
         Lymphatic_Endothelial_cell, Lymphoid, Mast_cell, Mesothelial_cell,
         Mural_cell, Myeloid, Neural_cell, Ventricular_Cardiomyocyte)
# Create a new column combining section and status
seurat_df <- seurat_df %>%
  mutate(section_status = paste0(section, “-”, status))
# Gather the data to long format
seurat_df_long <- seurat_df %>%
  gather(key = “cell_type”, value = “percentage”, -section, -status, -section_status)
# Calculate percentages
seurat_df_long <- seurat_df_long %>%
  group_by(section_status) %>%
  mutate(percentage = percentage / sum(percentage) * 100)
# Plot the stacked barplot with y-axis scaled to 0-100%
ggplot(seurat_df_long, aes(x = section_status, y = percentage, fill = cell_type)) +
  geom_bar(stat = “identity”, position = “fill”) +
  labs(title = “Percentage of Cell Types per Section and Status”,
       x = “Section and Status”,
       y = “Percentage”) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
###########################
# Figure 2B: Stacked plot of cell types
###########################
seurat_df <- spaceranger_integrated@meta.data %>%
  select(seurat_clusters, status, Adipocyte, Atrial_Cardiomyocyte, Endothelial_cell, Fibroblast,
         Lymphatic_Endothelial_cell, Lymphoid, Mast_cell, Mesothelial_cell,
         Mural_cell, Myeloid, Neural_cell, Ventricular_Cardiomyocyte)
# Create a new column combining seurat_clusters and status
seurat_df <- seurat_df %>%
  mutate(cluster_status = paste0(seurat_clusters, “-”, status))
# Gather the data to long format
seurat_df_long <- seurat_df %>%
  gather(key = “cell_type”, value = “percentage”, -seurat_clusters, -status, -cluster_status)
# Calculate percentages
seurat_df_long <- seurat_df_long %>%
  group_by(cluster_status) %>%
  mutate(percentage = percentage / sum(percentage) * 100)
# Plot the stacked barplot with y-axis scaled to 0-100%
ggplot(seurat_df_long, aes(x = cluster_status, y = percentage, fill = cell_type)) +
  geom_bar(stat = “identity”, position = “fill”) +
  labs(title = “Percentage of Cell Types per Seurat Cluster and Status”,
       x = “Seurat Cluster and Status”,
       y = “Percentage”) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
###########################
# Figure 2C: Stacked plot of cell types per cluster
###########################
seurat_df <- spaceranger_integrated@meta.data %>%
  select(seurat_clusters, Adipocyte, Atrial_Cardiomyocyte, Endothelial_cell, Fibroblast,
         Lymphatic_Endothelial_cell, Lymphoid, Mast_cell, Mesothelial_cell,
         Mural_cell, Myeloid, Neural_cell, Ventricular_Cardiomyocyte)
# Gather the data to long format
seurat_df_long <- seurat_df %>%
  gather(key = “cell_type”, value = “percentage”, -seurat_clusters)
# Calculate percentages
seurat_df_long <- seurat_df_long %>%
  group_by(seurat_clusters) %>%
  mutate(percentage = percentage / sum(percentage) * 100)
# Plot the stacked barplot with y-axis scaled to 0-100%
ggplot(seurat_df_long, aes(x = factor(seurat_clusters), y = percentage, fill = cell_type)) +
  geom_bar(stat = “identity”, position = “fill”) +
  labs(title = “Percentage of Cell Types per Seurat Cluster”,
       x = “Seurat Cluster”,
       y = “Percentage”) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0.35, linetype = “dashed”, color = “red”)
####################
# do the same excluding adipocytes and endothelial cells
seurat_df <- spaceranger_integrated@meta.data %>%
  select(seurat_clusters, Atrial_Cardiomyocyte, Fibroblast,
         Lymphatic_Endothelial_cell, Lymphoid, Mast_cell, Mesothelial_cell,
         Mural_cell, Myeloid, Neural_cell, Ventricular_Cardiomyocyte) %>%
  filter(seurat_clusters %in% c(3, 4, 6, 10, 11, 12, 13, 14))
# Gather the data to long format
seurat_df_long <- seurat_df %>%
  gather(key = “cell_type”, value = “percentage”, -seurat_clusters)
# Calculate percentages
seurat_df_long <- seurat_df_long %>%
  group_by(seurat_clusters) %>%
  mutate(percentage = percentage / sum(percentage) * 100)
# Plot the stacked barplot with y-axis scaled to 0-100%
ggplot(seurat_df_long, aes(x = factor(seurat_clusters), y = percentage, fill = cell_type)) +
  geom_bar(stat = “identity”, position = “fill”) +
  labs(title = “Percentage of Cell Types per Seurat Cluster”,
       x = “Seurat Cluster”,
       y = “Percentage”) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
######################
# make stacked bar plot sorting the y-axis based on adipocyte proportions
# Calculate the average proportion of Adipocytes for each cluster
adipocyte_proportion <- seurat_df_long %>%
  filter(cell_type == “Adipocyte”) %>%
  group_by(seurat_clusters) %>%
  summarise(adipocyte_proportion = mean(percentage)) %>%
  arrange(adipocyte_proportion)
# Reorder clusters based on the ascending average proportion of Adipocytes
seurat_df_long$seurat_clusters <- factor(seurat_df_long$seurat_clusters, levels = adipocyte_proportion$seurat_clusters)
# Plot the stacked barplot with y-axis scaled to 0-100%
ggplot(seurat_df_long, aes(x = seurat_clusters, y = percentage, fill = cell_type)) +
  geom_bar(stat = “identity”, position = “fill”) +
  labs(title = “Percentage of Cell Types per Seurat Cluster”,
       x = “Seurat Cluster”,
       y = “Percentage”) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
