###############################
# Run Pathway Finder
###############################
# lets try the pathway finder for DE genes from CSIDE
library(pathfindR)
library(biomaRt)

# first add gene symbols to the dataframes
# Specify the Ensembl dataset you want to use (human)
ensembl_dataset <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Define your Ensembl IDs
ensembl_ids_A1_adipo <- rownames(results_de_A1_adipo)
ensembl_ids_A1_endo <- rownames(results_de_A1_endo)
ensembl_ids_C1_adipo <- rownames(results_de_C1_adipo)
ensembl_ids_C1_endo <- rownames(results_de_C1_endo)

# Get gene symbols corresponding to Ensembl IDs
gene_symbols_A1_adipo <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                               filters = "ensembl_gene_id", 
                               values = ensembl_ids_A1_adipo, 
                               mart = ensembl_dataset)
gene_symbols_A1_endo <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                              filters = "ensembl_gene_id", 
                              values = ensembl_ids_A1_endo, 
                              mart = ensembl_dataset)
gene_symbols_C1_adipo <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                               filters = "ensembl_gene_id", 
                               values = ensembl_ids_C1_adipo, 
                               mart = ensembl_dataset)
gene_symbols_C1_endo <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                              filters = "ensembl_gene_id", 
                              values = ensembl_ids_C1_endo, 
                              mart = ensembl_dataset)

# Merge gene symbols with original data based on Ensembl IDs
merged_data_A1_adipo <- merge(gene_symbols_A1_adipo, results_de_A1_adipo, by.x = "ensembl_gene_id", by.y = "row.names")
merged_data_A1_endo <- merge(gene_symbols_A1_endo, results_de_A1_endo, by.x = "ensembl_gene_id", by.y = "row.names")
merged_data_C1_adipo <- merge(gene_symbols_C1_adipo, results_de_C1_adipo, by.x = "ensembl_gene_id", by.y = "row.names")
merged_data_C1_endo <- merge(gene_symbols_C1_endo, results_de_C1_endo, by.x = "ensembl_gene_id", by.y = "row.names")

# Extract relevant columns from merged data
merged_data_A1_adipo_gene_symbols <- merged_data_A1_adipo$external_gene_name
merged_data_A1_adipo_logFC <- merged_data_A1_adipo$`log_fc`
merged_data_A1_adipo_adj_P_Val <- merged_data_A1_adipo$`p_val`

merged_data_A1_endo_gene_symbols <- merged_data_A1_endo$external_gene_name
merged_data_A1_endo_logFC <- merged_data_A1_endo$`log_fc`
merged_data_A1_endo_adj_P_Val <- merged_data_A1_endo$`p_val`

merged_data_C1_adipo_gene_symbols <- merged_data_C1_adipo$external_gene_name
merged_data_C1_adipo_logFC <- merged_data_C1_adipo$`log_fc`
merged_data_C1_adipo_adj_P_Val <- merged_data_C1_adipo$`p_val`

merged_data_C1_endo_gene_symbols <- merged_data_C1_endo$external_gene_name
merged_data_C1_endo_logFC <- merged_data_C1_endo$`log_fc`
merged_data_C1_endo_adj_P_Val <- merged_data_C1_endo$`p_val`

# Create a dataframe with desired columns
result_df_A1_adipo <- data.frame(Gene.symbol = merged_data_A1_adipo_gene_symbols, logFC = merged_data_A1_adipo_logFC, adj.P.Val = merged_data_A1_adipo_adj_P_Val)
result_df_A1_endo <- data.frame(Gene.symbol = merged_data_A1_endo_gene_symbols, logFC = merged_data_A1_endo_logFC, adj.P.Val = merged_data_A1_endo_adj_P_Val)
result_df_C1_adipo <- data.frame(Gene.symbol = merged_data_C1_adipo_gene_symbols, logFC = merged_data_C1_adipo_logFC, adj.P.Val = merged_data_C1_adipo_adj_P_Val)
result_df_C1_endo <- data.frame(Gene.symbol = merged_data_C1_endo_gene_symbols, logFC = merged_data_C1_endo_logFC, adj.P.Val = merged_data_C1_endo_adj_P_Val)

###############################
# Do pathway analysis for up and down DE genes for adipocytes and endothelial cells of A1 and C1
###############################

output_pathways_A1_adipo_reac_up <- run_pathfindR(result_df_A1_adipo[result_df_A1_adipo$logFC > 0, ], gene_sets = "Reactome", p_val_threshold = 0.05)
output_pathways_A1_adipo_reac_down <- run_pathfindR(result_df_A1_adipo[result_df_A1_adipo$logFC < 0, ], gene_sets = "Reactome", p_val_threshold = 0.05)

output_pathways_A1_endo_reac_up <- run_pathfindR(result_df_A1_endo[result_df_A1_endo$logFC > 0, ], gene_sets = "Reactome", p_val_threshold = 0.05)
output_pathways_A1_endo_reac_down <- run_pathfindR(result_df_A1_endo[result_df_A1_endo$logFC < 0, ], gene_sets = "Reactome", p_val_threshold = 0.05)

output_pathways_C1_adipo_reac_up <- run_pathfindR(result_df_C1_adipo[result_df_C1_adipo$logFC > 0, ], gene_sets = "Reactome", p_val_threshold = 0.05)
output_pathways_C1_adipo_reac_down <- run_pathfindR(result_df_C1_adipo[result_df_C1_adipo$logFC < 0, ], gene_sets = "Reactome", p_val_threshold = 0.05)

output_pathways_C1_endo_reac_up <- run_pathfindR(result_df_C1_endo[result_df_C1_endo$logFC > 0, ], gene_sets = "Reactome", p_val_threshold = 0.05)
output_pathways_C1_endo_reac_down <- run_pathfindR(result_df_C1_endo[result_df_C1_endo$logFC < 0, ], gene_sets = "Reactome", p_val_threshold = 0.05)

# Save the object as an RDS file
output_path <- "/groups/umcg-franke-scrna/tmp01/projects/epifat/ongoing/seurat_preprocess_samples/objects/"
saveRDS(output_pathways_A1_adipo_reac_up, file.path(output_path, "output_pathways_A1_adipo_reac_up.rds"))
saveRDS(output_pathways_A1_adipo_reac_down, file.path(output_path, "output_pathways_A1_adipo_reac_down.rds"))
saveRDS(output_pathways_A1_endo_reac_up, file.path(output_path, "output_pathways_A1_endo_reac_up.rds"))
saveRDS(output_pathways_A1_endo_reac_down, file.path(output_path, "output_pathways_A1_endo_reac_down.rds"))
saveRDS(output_pathways_C1_adipo_reac_up, file.path(output_path, "output_pathways_C1_adipo_reac_up.rds"))
saveRDS(output_pathways_C1_adipo_reac_down, file.path(output_path, "output_pathways_C1_adipo_reac_down.rds"))
saveRDS(output_pathways_C1_endo_reac_up, file.path(output_path, "output_pathways_C1_endo_reac_up.rds"))
saveRDS(output_pathways_C1_endo_reac_down, file.path(output_path, "output_pathways_C1_endo_reac_down.rds"))

###############################
# load functions
###############################
split_string_with_max_length <- function(string, max_length=15) {
  # split the string on the space first
  words <- strsplit(string, ' ')[[1]]
  # and initial pasted string
  final_string <- ''
  # and the temporary string
  tmp_string <- ''
  # check each word
  for (word in words) {
    # get the lenght of the temporary string
    tmp_len <- nchar(tmp_string)
    # get the length of the word
    word_len <- nchar(word)
    # length of the final string
    final_len <- nchar(final_string)
    # check if adding the next word would make us go over the limit
    if ((tmp_len + word_len) > max_length) {
      if (final_len > 0) {
        # add the temporary word to the final string if it was not empty
        final_string <- paste(final_string, '\n', tmp_string, sep = '')
      }
      else{
        # if there was not final string yet, it will just become that temporary string
        final_string <- tmp_string
      }
      # reset the temporary string
      tmp_string <- word
    }
    else{
      # otherwise add it to the temporary string
      if (tmp_len > 0) {
        tmp_string <- paste(tmp_string, word)
      }
      # or set if the temporary string was empty
      else{
        tmp_string <- word
      }
    }
  }
  # add last word
  if (nchar(paste(final_string, tmp_string)) > max_length) {
    if (nchar(final_string) == 0) {
      final_string <- tmp_string
    }
    else if (nchar(tmp_string > 0)) {
      final_string <- paste(final_string, '\n', tmp_string, sep = '')
    }
    else{
      # nothing the final string is as it is
    }
  }
  else {
    if (nchar(final_string) == 0) {
      final_string <- tmp_string
    }
    else if(nchar(tmp_string) == 0){
      # nothing, the final string is as is
    }
    else{
      final_string <- paste(final_string, ' ', tmp_string, sep = '')
    }
  }
  # add the last word
  return(final_string)
}

# Get Reactome pathways for Treeplot
get_children <- function(relation_table, starting_id){
  # get all of the children of the starting ID
  children <- as.character(relation_table[relation_table$V1 == starting_id, 'V2'])
  # these children are all family
  family <- children
  # see if there were any children
  if(length(children) > 0){
    # if there were children, we need to get their children as well
    for(child in children){
      # get the grandchildren and add these to the family
      grand_children <- get_children(relation_table, child)
      family <- c(family, grand_children)
    }
  }
  return(family)
}

get_filtered_pathway_names <- function(pathway_table, relation_table, starting_id, add_id=F, only_id=F){
  # get all of the children of the starting ID
  all_children <- get_children(relation_table, starting_id)
  # get the names of the pathways that are children
  pathways_to_include_names <- pathway_table[pathway_table$V1 %in% all_children, ]
  pathway_names <- NULL
  # add the id if requested
  if (only_id) {
    pathway_names <-pathways_to_include_names[['V1']]
  }
  else if (add_id) {
    pathway_names <- paste(pathways_to_include_names[['V2']], pathways_to_include_names[['V1']])
  }
  else {
    pathway_names <-pathways_to_include_names[['V2']]
  }
  return(pathway_names)
}

pathway_names_to_top_level <- function(pathway_mappings, pathway_names) {
  # subset to just human to speed up the search
  pathway_names <- pathway_names[pathway_names$V3 == 'Homo sapiens', ]
  # get the top level human pathways
  top_level_pathways <- unique(pathway_mappings[pathway_mappings$V2 %in% pathways$V1 & # is human
                                                  !(pathway_mappings$V1 %in% pathway_mappings$V2) , 'V1']) # doesn't have a parent
  # first store in list
  pathway_name_to_top_level <- list()
  # check each top level pathway
  for (top_level in top_level_pathways) {
    # get the filtered names for that one
    filtered_names <- get_filtered_pathway_names(pathways, pathway_mappings, top_level, add_id = T)
    # and filtered IDs
    filtered_ids <- get_filtered_pathway_names(pathways, pathway_mappings, top_level, only_id = T)
    # get the name of this parent
    parent_name <- pathway_names[match(top_level, pathway_names[['V1']]), 'V2']
    # make dataframe
    pathway_mapping_category <- data.frame(pathway = filtered_names, id = filtered_ids, top_parent = rep(parent_name, times = length(parent_name)), top_id = rep(top_level, length(parent_name)))
    # add the category itself as well
    pathway_mapping_category <- rbind(pathway_mapping_category, data.frame(pathway = c(paste(parent_name, top_level)), id = c(top_level), top_parent = c(parent_name), top_id=c(top_level)))
    # add to list
    pathway_name_to_top_level[[top_level]] <- pathway_mapping_category
  }
  # merge lists
  pathways_to_top_parent <- do.call('rbind', pathway_name_to_top_level)
  return(pathways_to_top_parent)
}

# treeplot function
plot_pathway_categories_celltype <- function(pathway_table, category_column='category', legendless=F, paper_style=T, to_fractions=F, use_distinct_colours=T, use_sampling=F, color_indices=NULL, drop_categories=NULL, text_loc='topleft') {
  # now turn this into a table
  number_per_category <- data.frame(table(pathway_table_significant[[category_column]]))
  # set better column names
  colnames(number_per_category) <- c('category', 'number')
  # turn into plot
  p <- ggplot(number_per_category, aes(area = number, fill = category, label = category)) +
    geom_treemap() +
    geom_treemap_text(place = text_loc)
  # use distinct colours if requested
  if (use_distinct_colours) {
    # get the unique possible assigments
    possible_assignments <- unique(pathway_table_significant[['category']])
    # get an equal amount of colours
    possible_colours <- NULL
    if(length(possible_assignments) > 74) {
      possible_colours <- sample_tons_of_colors(length(possible_assignments), use_sampling = use_sampling, color_indices = color_indices)
    }
    else {
      possible_colours <- sample_many_colours(length(possible_assignments), use_sampling = use_sampling, color_indices = color_indices)
    }
    # put into a list
    colour_mapping <- as.list(possible_colours)
    names(colour_mapping) <- possible_assignments
    # add to plot
    p <- p + scale_fill_manual(values = colour_mapping)
  }
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  if (paper_style) {
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  return(p)
}

# treeplot colours
#' get a vector of as distinct possible colours
#' 
#' @param number_of_colours how many colours to return
#' @param use_sampling whether or not to randomly extract the colours instead of grabbing the first n colours
#' @param color_indices (optional, not used by default) if specific colours are needed, supply the indices of the colours here. Use 'get_available_colours_grid' to get the colours and their indices
#' @returns a vector of colours
#' 
sample_many_colours <- function(number_of_colours, use_sampling=F, color_indices=NULL) {
  # get all colours from the 'quality' palettes
  quality_colour_palettes <- brewer.pal.info[brewer.pal.info[['category']] == 'qual', ]
  # save each palette
  colours_per_palette <- list()
  # apply over each palette
  for (i in 1:nrow(quality_colour_palettes)) {
    # get the name of the palette
    palette_name <- rownames(quality_colour_palettes)[i]
    # get the number of colours in the palette
    palette_max_colours <- quality_colour_palettes[i, 'maxcolors']
    # use brewer.pal to get all colours
    colours_palette <- brewer.pal(palette_max_colours, palette_name)
    # put result in the list
    colours_per_palette[[palette_name]] <- colours_palette
  }
  # merge all palettes
  all_colours <- do.call('c', colours_per_palette)
  # randomly get colours from the palette
  max_possible_colours <- length(all_colours)
  if (is.null(number_of_colours)) {
    message('no number of colors supplied, assuming color indices have been')
  }
  else if (number_of_colours > max_possible_colours) {
    message(paste('requesting more colours than is possible: ', as.character(number_of_colours), ' vs ', max_possible_colours, ', returning max possible', sep = ''))
    number_of_colours <- max_possible_colours
  }
  colours_to_return <- NULL
  # specific colours we like (the indices)
  if (!is.null(color_indices)) {
    colours_to_return <- all_colours[color_indices]
  }
  # or use sampling
  else if (use_sampling) {
    colours_to_return <- sample(all_colours, number_of_colours)
  }
  # or the first x colours
  else {
    colours_to_return <- all_colours[1 : number_of_colours]
  }
  return(colours_to_return)
}

#' get a vector of as distinct possible colours, but with more possibilities ()
#' 
#' @param number_of_colours how many colours to return
#' @param use_sampling whether or not to randomly extract the colours instead of grabbing the first n colours
#' @param color_indices (optional, not used by default) if specific colours are needed, supply the indices of the colours here. Use 'get_available_colours_grid' to get the colours and their indices
#' @returns a vector of colours
#' 
sample_tons_of_colors <- function(number_of_colours, use_sampling=F, color_indices=NULL) {
  # get colours available to device
  all_colours <- grDevices::colors()
  # remove gray
  all_colours <- all_colours[grep('gr(a|e)y', all_colours, invert = T)]
  # check how many are possible
  max_possible_colours <- length(all_colours)
  if (is.null(number_of_colours)) {
    message('no number of colors supplied, assuming color indices have been')
  }
  else if (number_of_colours > max_possible_colours) {
    message(paste('requesting more colours than is possible: ', as.character(number_of_colours), ' vs ', max_possible_colours, ', returning max possible', sep = ''))
    number_of_colours <- max_possible_colours
  }
  colours_to_return <- NULL
  # specific colours we like (the indices)
  if (!is.null(color_indices)) {
    colours_to_return <- all_colours[color_indices]
  }
  # or use sampling
  else if (use_sampling) {
    colours_to_return <- sample(all_colours, number_of_colours)
  }
  # or the first x colours
  else {
    colours_to_return <- all_colours[1 : number_of_colours]
  }
  return(colours_to_return)
}


#' get a grid showing the available colours and their indices
#' 
#' @param many use the 'many' method to get the colours
#' @param tons use the 'tons' method to get the colours
#' @returns a ggplot grid showing the available colours and their indices
#' 
get_available_colours_grid <- function(many=T, tons=F) {
  colours_possible <- NULL
  # get from the many method
  if (many) {
    # ask for unreasonable amount
    colours_possible <- sample_many_colours(1000)
  }
  else if(tons) {
    colours_possible <- sample_tons_of_colors(1000)
  }
  # get how many colours we actually have
  available_colours <- length(colours_possible)
  # we need to put that into a square grid, so we need to get the square root, to know how many rows and columns
  nrow_and_ncol <- sqrt(available_colours)
  # and we need to round that up of course
  nrow_and_ncol <- ceiling(nrow_and_ncol)
  # so we'll have a total number of blocks
  total_cells <- nrow_and_ncol * nrow_and_ncol
  # let's see how many colours we are off from that number of cells
  cells_no_colour_number <- total_cells - available_colours
  # we will just add white for those
  cells_no_colour <- rep('white', times = cells_no_colour_number)
  # add that to the colours we have
  colours_possible <- c(colours_possible, cells_no_colour)
  # create each combination of x and y
  indices_grid <- expand.grid(as.character(1 : nrow_and_ncol), as.character(1 : nrow_and_ncol))
  # add the index and colour name
  indices_grid[['index_colour']] <- paste(c(1:total_cells), colours_possible, sep = '\n')
  # make mapping of colours
  colours_to_use <- as.list(colours_possible)
  names(colours_to_use) <- indices_grid[['index_colour']]
  # now plot
  p <- ggplot(data = indices_grid, mapping = aes(x = Var1, y = Var2, fill = index_colour)) + 
    geom_tile() + 
    geom_text(aes(label=index_colour)) + 
    scale_fill_manual(values = colours_to_use) + 
    theme(legend.position = 'none')
  return(p)
}

plot_pathway_categories_celltype <- function(pathway_table_significant, category_column='category', legendless=F, paper_style=T, to_fractions=F, use_distinct_colours=T, use_sampling=F, color_indices=NULL, drop_categories=NULL, text_loc='topleft') {
  # set NA to unannotated
  pathway_table_significant[is.na(pathway_table_significant[[category_column]]), category_column] <- 'unannotated'
  # remove categories if requested
  if (!is.null(drop_categories)) {
    pathway_table_significant <- pathway_table_significant[!(pathway_table_significant[[category_column]] %in% drop_categories), ]
  }
  # now turn this into a table
  number_per_category <- data.frame(table(pathway_table_significant[[category_column]]))
  # set better column names
  colnames(number_per_category) <- c('category', 'number')
  # turn into plot
  p <- ggplot(number_per_category, aes(area = number, fill = category, label = category)) +
    geom_treemap() +
    geom_treemap_text(place = text_loc)
  # use distinct colours if requested
  if (use_distinct_colours) {
    # get the unique possible assigments
    possible_assignments <- unique(pathway_table_significant[['category']])
    # get an equal amount of colours
    possible_colours <- NULL
    if(length(possible_assignments) > 74) {
      possible_colours <- sample_tons_of_colors(length(possible_assignments), use_sampling = use_sampling, color_indices = color_indices)
    }
    else {
      possible_colours <- sample_many_colours(length(possible_assignments), use_sampling = use_sampling, color_indices = color_indices)
    }
    # put into a list
    colour_mapping <- as.list(possible_colours)
    names(colour_mapping) <- possible_assignments
    # add to plot
    p <- p + scale_fill_manual(values = colour_mapping)
  }
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  if (paper_style) {
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  return(p)
}

plot_pathway_categories_celltype <- function(pathway_table_significant, category_column='category', legendless=F, paper_style=T, to_fractions=F, use_distinct_colours=T, use_sampling=F, color_indices=NULL, drop_categories=NULL, text_loc='topleft') {
  # set NA to unannotated
  pathway_table_significant[is.na(pathway_table_significant[[category_column]]), category_column] <- 'unannotated'
  # remove categories if requested
  if (!is.null(drop_categories)) {
    pathway_table_significant <- pathway_table_significant[!(pathway_table_significant[[category_column]] %in% drop_categories), ]
  }
  # now turn this into a table
  number_per_category <- data.frame(table(pathway_table_significant[[category_column]]))
  # set better column names
  colnames(number_per_category) <- c('category', 'number')
  # turn into plot
  p <- ggplot(number_per_category, aes(area = number, fill = category, label = category)) +
    geom_treemap() +
    geom_treemap_text(place = text_loc)
  # use distinct colours if requested
  if (use_distinct_colours) {
    # get the unique possible assigments
    possible_assignments <- unique(number_per_category[['category']])
    # get an equal amount of colours
    possible_colours <- NULL
    if(length(possible_assignments) > 74) {
      possible_colours <- sample_tons_of_colors(length(possible_assignments), use_sampling = use_sampling, color_indices = color_indices)
    }
    else {
      possible_colours <- sample_many_colours(length(possible_assignments), use_sampling = use_sampling, color_indices = color_indices)
    }
    # put into a list
    colour_mapping <- as.list(possible_colours)
    names(colour_mapping) <- possible_assignments
    # add to plot
    p <- p + scale_fill_manual(values = colour_mapping)
  }
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  if (paper_style) {
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  return(p)
}

#########################
# load libraries
#########################
library(RColorBrewer)
library(treemapify) # not in container

#########################
# do stuff
#########################
# load the pathways
pathways <- read.table('/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/pathway_enrichment/ReactomePathways.tsv', sep='\t', quote = '')
# subset to just human to speed up the search
pathways <- pathways[pathways$V3 == 'Homo sapiens', ]
# load the pathway mapping
pathway_mappings <- read.table('/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/pathway_enrichment/ReactomePathwaysRelation.tsv', sep = '\t')
# get the filtered names
filtered_names <- c(get_filtered_pathway_names(pathways, pathway_mappings, 'R-HSA-168256', add_id = T), get_filtered_pathway_names(pathways, pathway_mappings, 'R-HSA-162582', add_id = T))
filtered_ids <- c(get_filtered_pathway_names(pathways, pathway_mappings, 'R-HSA-168256', only_id = T), get_filtered_pathway_names(pathways, pathway_mappings, 'R-HSA-162582', only_id = T))
# get the top level human pathways
pathways_to_top_level <- pathway_names_to_top_level(pathway_mappings = pathway_mappings, pathway_names = pathways)
pathways_to_top_level[is.na(pathways_to_top_level[['top_parent']]), 'top_parent'] <- 'unannotated'
# make a copy with more readable names
pathways_to_top_level[['top_parent_readable']] <- apply(pathways_to_top_level, 1, function(x){split_string_with_max_length(x[['top_parent']])})

# match ID results to id for up and down DE genes
output_pathways_A1_adipo_reac_up[["catagory"]] <- pathways_to_top_level[match(output_pathways_A1_adipo_reac_up[["ID"]], pathways_to_top_level[["id"]]), "top_parent"]
output_pathways_A1_adipo_reac_down[["catagory"]] <- pathways_to_top_level[match(output_pathways_A1_adipo_reac_down[["ID"]], pathways_to_top_level[["id"]]), "top_parent"]
output_pathways_A1_endo_reac_up[["catagory"]] <- pathways_to_top_level[match(output_pathways_A1_endo_reac_up[["ID"]], pathways_to_top_level[["id"]]), "top_parent"]
output_pathways_A1_endo_reac_down[["catagory"]] <- pathways_to_top_level[match(output_pathways_A1_endo_reac_down[["ID"]], pathways_to_top_level[["id"]]), "top_parent"]
output_pathways_C1_adipo_reac_up[["catagory"]] <- pathways_to_top_level[match(output_pathways_C1_adipo_reac_up[["ID"]], pathways_to_top_level[["id"]]), "top_parent"]
output_pathways_C1_adipo_reac_down[["catagory"]] <- pathways_to_top_level[match(output_pathways_C1_adipo_reac_down[["ID"]], pathways_to_top_level[["id"]]), "top_parent"]
output_pathways_C1_endo_reac_up[["catagory"]] <- pathways_to_top_level[match(output_pathways_C1_endo_reac_up[["ID"]], pathways_to_top_level[["id"]]), "top_parent"]
output_pathways_C1_endo_reac_down[["catagory"]] <- pathways_to_top_level[match(output_pathways_C1_endo_reac_down[["ID"]], pathways_to_top_level[["id"]]), "top_parent"]

# look at colors
get_available_colours_grid()

# make treeplot for up and down DE genes
A1_adipo_up <- plot_pathway_categories_celltype(output_pathways_A1_adipo_reac_up, category_column='catagory', use_distinct_colours=T)
A1_adipo_down <- plot_pathway_categories_celltype(output_pathways_A1_adipo_reac_down, category_column='catagory', use_distinct_colours=T)
A1_endo_up <- plot_pathway_categories_celltype(output_pathways_A1_endo_reac_up, category_column='catagory', use_distinct_colours=T)
A1_endo_down <- plot_pathway_categories_celltype(output_pathways_A1_endo_reac_down, category_column='catagory', use_distinct_colours=T)
C1_adipo_up <- plot_pathway_categories_celltype(output_pathways_C1_adipo_reac_up, category_column='catagory', use_distinct_colours=T)
C1_adipo_down <- plot_pathway_categories_celltype(output_pathways_C1_adipo_reac_down, category_column='catagory', use_distinct_colours=T)
C1_endo_up <- plot_pathway_categories_celltype(output_pathways_C1_endo_reac_up, category_column='catagory', use_distinct_colours=T)
C1_endo_down <- plot_pathway_categories_celltype(output_pathways_C1_endo_reac_down, category_column='catagory', use_distinct_colours=T)
# adipo up --> signal transduction
# adipo down --> metabolism
# endo up --> immune system
# endo down --> immune system

# filter in up and down pathways
include_pathways <- pathways_to_top_level[pathways_to_top_level[['top_parent']] == 'Immune system', 'id']
new_A1_adipo_up <- output_pathways_A1_adipo_reac_up[output_pathways_A1_adipo_reac_up[['catagory']] %in% c('Signal Transduction'), ] 
enrichment_chart(new_A1_adipo_up, plot_by_cluster = T)
new_A1_adipo_down <- output_pathways_A1_adipo_reac_down[output_pathways_A1_adipo_reac_down[['catagory']] %in% c('Metabolism'), ] 
enrichment_chart(new_A1_adipo_down, plot_by_cluster = T)
new_A1_endo_up <- output_pathways_A1_endo_reac_up[output_pathways_A1_endo_reac_up[['catagory']] %in% c('Immune System'), ] 
enrichment_chart(new_A1_endo_up, plot_by_cluster = T)
new_A1_endo_down <- output_pathways_A1_endo_reac_down[output_pathways_A1_endo_reac_down[['catagory']] %in% c('Immune System'), ] 
enrichment_chart(new_A1_endo_down, plot_by_cluster = T)
new_C1_adipo_up <- output_pathways_C1_adipo_reac_up[output_pathways_C1_adipo_reac_up[['catagory']] %in% c('Signal Transduction'), ] 
enrichment_chart(new_C1_adipo_up, plot_by_cluster = T)
new_C1_adipo_down <- output_pathways_C1_adipo_reac_down[output_pathways_C1_adipo_reac_down[['catagory']] %in% c('Metabolism'), ] 
enrichment_chart(new_C1_adipo_down, plot_by_cluster = T)
new_C1_endo_up <- output_pathways_C1_endo_reac_up[output_pathways_C1_endo_reac_up[['catagory']] %in% c('Signal Transduction'), ] 
enrichment_chart(new_C1_endo_up, plot_by_cluster = T)
new_C1_endo_down <- output_pathways_C1_endo_reac_down[output_pathways_C1_endo_reac_down[['catagory']] %in% c('Immune System'), ] 
enrichment_chart(new_C1_endo_down, plot_by_cluster = T)



