

library(tidyverse)
# library(pbapply)
# library(readxl)
# 
# library(tidygraph)
# library(ggraph)
# library(ggrepel)
# library(ggthemes)
# library(ggrastr)
# library(ggforce)
# library(ggupset)
# library(plotly)
# library(patchwork)
# 
# library(uwot)

source("scripts/functions_utility.R")

select <- dplyr::select

dataset_metadata <- read_csv("run_settings/20210928_settings.csv")

# data_paths <- read_csv("run_settings/20211012 all_datasets.csv")

file_structure <-
  create_folder_structure(dataset_id = dataset_metadata$dataset_id,
                          dataset_run_id = dataset_metadata$version,
                          main_folder = "results/Clustering_results", # CHANGE TO AVOID SPACES
                          folders = c("svg" = "svg",
                                      "data" = "data",
                                      "treemap" = "svg/treemap",
                                      "heatmap" = "svg/heatmap",
                                      "PCA" = "PCA",
                                      "distance" = "distance",
                                      "graph" = "graph",
                                      "UMAP" = "UMAP",
                                      "enrichment" = "enrichment",
                                      "clustering" = "clustering",
                                      "evaluation" = "evaluation"))


a <- 
  readRDS("results/Clustering_results/singlecell_HPA21v1/graph/neighbors.rds")
a_dist <- 
  readRDS("results/Clustering_results/singlecell_HPA21v1/distance/distances.rds")

a_dist_m <- 
  as.matrix(a_dist)

ex_gene <- 
  colnames(a$nn)[1]

other_genes1 <- 
  a$nn[ex_gene, ] %>% 
  {.[which(. == 1)]} %>% 
  names()

other_genes2 <- 
  a$nn[, ex_gene] %>% 
  {.[which(. == 1)]} %>% 
  names()


gene_nn %>%
  mutate(gen1 = ifelse(gene_id < name, 
                       gene_id,
                       name),
         gen2 = ifelse(gene_id < name, 
                       name,
                       gene_id)) %>% 
  select(gen1, gen2, value) %>% 
  distinct() %>% 
  ggplot(aes(value)) +
  geom_density()

distances <- 
  readRDS(file_) %>% 
  as.matrix()

nn_data_dist <- 
  colnames(a_dist_m) %>% 
  enframe("i", "gene") %>% 
  select(gene) %>% 
  group_by(gene) %>% 
  do({
    g_gene <<- .$gene 
    a_dist_m[g_gene, ] %>% 
      enframe("neighbor_gene", "cor_distance") %>% 
      arrange(cor_distance) %>% 
      head(20)
  }) %>% 
  ungroup()

nn_data_nn <- 
  colnames(a$nn) %>% 
  enframe("i", "gene") %>% 
  select(gene) %>% 
  group_by(gene) %>% 
  do({
    g_gene <<- .$gene 
    a$nn[g_gene, ] %>% 
      {.[which(. == 1)]} %>% 
      names() %>% 
      enframe
  }) %>% 
  ungroup()


nn_data_dist$neighbor_gene == nn_data_nn$value
nn_data_dist$neighbor_gene == nn_data_nn$value
nn_data_dist %>% 
  inner_join(nn_data_nn,
             by = c("gene", "neighbor_gene" = "value"))  
  filter(neigh)
####



neighbor_data <- 
  file_structure %>% 
  map(function(x) x$distance) %>% 
  map(function(x) {
    file_ <- paste(x, "distances.rds", sep = "/")
    if(file.exists(file_)) {
      distances <- 
        readRDS(file_) %>% 
        as.matrix()
      
      nn_data <- 
        colnames(distances) %>% 
        enframe("i", "gene") %>% 
        select(gene) %>% 
        group_by(gene) %>% 
        do({
          g_gene <<- .$gene 
          distances[g_gene, ] %>% 
            enframe("neighbor_gene", "cor_distance") %>% 
            arrange(cor_distance) %>% 
            head(20)
        }) %>% 
        ungroup()
      rm(distances)
      
      nn_data
      
    } else {
      NULL
    }
  }) %>% 
  bind_rows(.id = "dataset_id") 



gene_nn <-
  colnames(a_dist_m) %>% 
  enframe("i", "gene_id") %>% 
  select(gene_id) %>% 
  group_by(gene_id) %>% 
  do({
    g_gene <- .$gene_id 
    a_dist_m[g_gene, ] %>% 
      enframe("neighbor_gene", "cor_distance") %>% 
      arrange(value) %>% 
      filter(name != g_gene) %>% 
      head(20)
  }) %>% 
  ungroup()



