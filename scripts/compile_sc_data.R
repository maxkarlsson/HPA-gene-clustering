

library(tidyverse)


sc_data <- 
  read_tsv("../../Data/HPA/HPA22_E103_files/single cell/single_cell_cluster_data_v3_103.tsv")

sc_anno <- 
  read_tsv("../../Data/HPA/HPA22_E103_files/single cell/single_cell_annotation_v3_103.tsv")


sc_data_formatted <-
  sc_data %>% 
  left_join(sc_anno) %>% 
  select(ensg_id, assay_id, tissue_name, cell_type_name, cluster_id, ntpm) %>% 
  mutate(sample_identifyer = paste(assay_id, cluster_id, sep = "_")) %>% 
  unite(sample, tissue_name, cell_type_name, cluster_id) %>% 
  unite(sample, sample_identifyer, sample, sep = ";") %>% 
  select(ensg_id, sample, ntpm)

sc_data_formatted %>% 
  write_tsv("data/expression_data/HPA22_E103/sc_cluster_data.tsv")


