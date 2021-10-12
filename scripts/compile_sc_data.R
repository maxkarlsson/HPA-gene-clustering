

library(tidyverse)


sc_data <- 
  read_tsv("../../Data/HPA/HPA21_E103_files/single cell/cluster_data.tsv")

sc_anno <- 
  read_tsv("../../Data/HPA/HPA21_E103_files/single cell/annotation.tsv")


sc_data_formatted <- 
  sc_data %>% 
  left_join(sc_anno) %>% 
  select(ensg_id, tissue_name, cell_type_name, cluster_id, ntpm) %>% 
  unite(sample, tissue_name, cell_type_name, cluster_id)

sc_data_formatted %>% 
  write_tsv("data/expression_data/HPA21_E103/sc_cluster_data.tsv")


