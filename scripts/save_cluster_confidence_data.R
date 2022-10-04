
# Save cluster confidence data from tissue expression clusters (HPAv21)

library(tidyverse)

source("scripts/functions_utility.R")

dataset_metadata <- read_csv("run_settings/20210928_settings.csv") #HPA21 v1

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

# Tissue memebership
membership_file <- 
  paste(file_structure[["tissue"]]$clustering, "cluster_memberships.tsv", sep = "/")

cluster_confidence_tissue <- 
  read_tsv(membership_file) %>% 
  as_tibble() %>% 
  rename(cluster_confidence =  membership) 

cluster_confidence_tissue_wide <- 
  cluster_confidence_tissue %>% 
  spread(cluster,cluster_confidence, fill = 0, sep = "_")

saveRDS(cluster_confidence_tissue_wide, "data/processed/cluster_confidence_tissue.rds")
write_csv(cluster_confidence_tissue_wide, "data/processed/cluster_confidence_tissue.csv")


# plot
library(pheatmap)
cluster_confidence_tissue %>% 
  spread(cluster,cluster_confidence, fill = 0, sep = "_") %>% 
  column_to_rownames("gene") %>% 
  pheatmap()

#Check

readRDS("data/processed/cluster_confidence_tissue.rds")
a <- read_csv("data/processed/cluster_confidence_tissue.csv")

a %>% gather(cluster, conf, -gene) %>% group_by(gene) %>% summarise(s_c = sum(conf)) %>% arrange(s_c)

#Single cell membership
membership_file <- 
  paste(file_structure[["singlecell"]]$clustering, "cluster_memberships.tsv", sep = "/")

cluster_confidence_singlecell <- 
  read_tsv(membership_file) %>% 
  as_tibble() %>% 
  rename(cluster_confidence =  membership) 

cluster_confidence_singlecell_wide <- 
  cluster_confidence_singlecell %>% 
  spread(cluster,cluster_confidence, fill = 0, sep = "_")

saveRDS(cluster_confidence_singlecell_wide, "data/processed/cluster_confidence_singlecell.rds")
write_csv(cluster_confidence_singlecell_wide, "data/processed/cluster_confidence_singlecell.csv")


# plot
library(pheatmap)
cluster_confidence_singlecell %>% 
  spread(cluster,cluster_confidence, fill = 0, sep = "_") %>% 
  column_to_rownames("gene") %>% 
  pheatmap()

#Check

readRDS("data/processed/cluster_confidence_singlecell.rds")
a <- read_csv("data/processed/cluster_confidence_singlecell.csv")

a %>% gather(cluster, conf, -gene) %>% group_by(gene) %>% summarise(s_c = sum(conf)) %>% arrange(-s_c)


# Some genes
a %>% 
  filter(gene == "ENSG00000176105") %>%
  gather(cluster, conf, -1) %>% 
  arrange(-conf)
