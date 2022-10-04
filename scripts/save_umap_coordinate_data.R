
# Save UMAP coordinates from single cell expression clusters 

library(tidyverse)

source("scripts/functions_utility.R")

dataset_metadata <- read_csv("run_settings/20220222_settings.csv")

data_paths <- read_csv("run_settings/20211012 all_datasets.csv")

file_structure <-
  create_folder_structure(dataset_id = dataset_metadata$dataset_id,
                          dataset_run_id = dataset_metadata$version,
                          main_folder = "results/Clustering_results",
                          folders = c("data" = "data", 
                                      "UMAP" = "UMAP",
                                      "svg" = "svg",
                                      "bubbleheatmap" = "svg/bubble",
                                      "treemap" = "svg/treemap",
                                      "heatmap" = "svg/heatmap",
                                      "enrichment" = "enrichment",
                                      "clustering" = "clustering"))
umap_data <- 
  file_structure %>% 
  map(function(x) x$UMAP) %>% 
  map(function(x) {
    file_ <- paste(x, "UMAP.tsv", sep = "/")
    if(file.exists(file_)) {
      read_tsv(file_)
    } else {
      NULL
    }
  }) %>% 
  bind_rows(.id = "dataset_id") 

clustering_data <- 
  file_structure %>% 
  map(function(x) x$clustering) %>% 
  map(function(x) {
    file_ <- paste(x, "final_consensus.tsv", sep = "/")
    if(file.exists(file_)) {
      read_tsv(file_)
    } else {
      NULL
    }
  }) %>% 
  bind_rows(.id = "dataset_id") %>% 
  mutate(cluster = as.character(cluster)) 


sc_gene_umap <- 
  umap_data %>% 
  filter(dataset_id == "singlecell") %>% 
  select(-dataset_id, -UMAP_1_scaled, -UMAP_2_scaled)

sc_gene_clusters_umap <- 
  umap_data %>% 
  filter(dataset_id == "singlecell") %>% 
  select(-dataset_id, -UMAP_1_scaled, -UMAP_2_scaled) %>% 
  left_join(clustering_data %>% 
              filter(dataset_id =="singlecell") %>% 
              select(-dataset_id) %>% 
              rename(expression_cluster = cluster))

saveRDS(sc_gene_umap, "data/processed/sc_gene_umap.rds")
saveRDS(sc_gene_clusters_umap, "data/processed/sc_gene_clusters_umap.rds")


sc_gene_clusters_umap %>% 
  ggplot(aes(UMAP_1,UMAP_2,color =expression_cluster)) +
  geom_point(show.legend = F) +
  theme_void()


readRDS("data/processed/sc_gene_clusters_umap.rds")

