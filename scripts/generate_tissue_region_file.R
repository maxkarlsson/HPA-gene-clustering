


library(tidyverse)


brain_region <- 
  read_tsv("../../Data/HPA/HPA21_E103_files/For gene clustering/consensus_brain_103.tsv")


tissue_region <- 
  read_tsv("../../Data/HPA/HPA21_E103_files/For gene clustering/consensus_aggregated_tissues_103.tsv")

tissue_hierarchy <- 
  read_tsv("../../Data/HPA/HPA21_E103_files/final_consensus_hierarchy_103.tsv")

brain_regions <- 
  tissue_hierarchy %>% 
  filter(organ == "Brain") %>% 
  pull(tissue)

tissue_region %>% 
  select(from_tissue) %>% 
  distinct()

tissue_region %>% 
  select(1, tissue = 2, ntpm) %>%
  filter(!tissue %in% brain_regions) %>% 
  bind_rows(brain_region %>% 
              select(ensg_id, tissue, ntpm)) %>% 
  write_tsv("../../Data/HPA/HPA21_E103_files/For gene clustering/region_aggregated_tissue_103.tsv")
  
