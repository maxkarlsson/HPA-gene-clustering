

library(tidyverse)
library(vroom)
library(magrittr)
source("scripts/functions_utility.R")



meta_data <- 
  read_csv("../../Data/HPA/HPA_metadata.csv")

meta_data_sc <- 
  read_tsv("data/meta/cell_type_annotation_200628_complete.tab") 

# Read example data if file exists. Otherwise, compile the data file (works only on Max's computer)

save_id <- 
  "all_cluster_data_mockup2.rds"
save_id_meta <- 
  "all_cluster_data_metadata_mockup2.rds"


all_data <- 
  compile_example_data(save_id)

dataset_metadata <-
  all_data %>% 
  names() %>% 
  enframe("i", "dataset_id") %>% 
  mutate(dataset = gsub("_.*", "", dataset_id),
         type = case_when(dataset_id %in% c("tissue", "blood_consensus", "celline_consensus") ~ "region",
                          dataset_id %in% c("brain") ~ "tissue",
                       T ~ "sample")) %>% 
  select(-i)

  
singlecell_metadata <-
  vroom("data/meta/rna_single_cell_type_tissue.tsv") %>%
  select(tissue_name = Tissue, 
         cluster = Cluster,
         cell_type = `Cell type`) %>% 
  distinct() %>% 
  unite(id, tissue_name, cluster, remove = F)

all_data_meta <-
  all_data %>% 
  map( . %>% 
         select(-1) %>% 
         names() %>% 
         enframe("i", "sample_id")) %>% 
  bind_rows(.id = "dataset_id") %>% 
  select(-i) %>%
  left_join(dataset_metadata) %>%
  left_join(meta_data %>% 
              select(sample, scilifelab_id),
            by = c("sample_id" = "scilifelab_id")) %>% 
  # filter(type != "consensus") %>%  
  left_join(singlecell_metadata, 
            by = c("sample_id" = "id")) %>% 
  mutate(sample = case_when(dataset == "singlecell" ~ cell_type,
                            type == "sample" ~ sample, 
                            type == "tissue" ~ sample_id, 
                            type == "region" ~ sample_id)) %>% 
  select(dataset_id, dataset, type, sample_id, sample)

saveRDS(list(all_data_meta = all_data_meta,
             dataset_metadata = dataset_metadata), 
        paste0("data/processed/", save_id_meta))
