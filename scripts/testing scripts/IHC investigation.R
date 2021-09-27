

IHC_annotation <- 
  read_tsv("data/meta/IHC_annotation_normal_tissue.tsv")


IHC_annotation %>% 
  filter(Level != "Not detected") %>% 
  filter(Reliability != "Uncertain")  %>% 
  select(Tissue, `Cell type`) %>% 
  distinct() %>% 
  group_by(`Cell type`) %>% 
  summarise(n = n_distinct(Tissue), 
            Tissue = paste(Tissue, collapse = ";")) %>% 
  arrange(-n) %>% 
  filter(n > 1) %>% 
  write_csv("results/IHC doublet cell names.csv")


IHC_annotation %>% 
  filter(Level %in% c("Low", "Medium", "High")) %>% 
  filter(Reliability != "Uncertain") %>% 
  group_by(Tissue, `Cell type`) %>% 
  count %>% 
  arrange(-n) %>% 
  print(n = 100)


IHC_annotation %>% 
  filter(Level %in% c("Low", "Medium", "High")) %>% 
  filter(Reliability != "Uncertain") %>% 
  select(Gene, `Cell type`) %>% 
  group_by(Gene) %>% 
  count() %>% 
  group_by(n) %>% 
  count %>% 
  filter(n <= 10) %>% 
  pull(nn) %>% 
  sum

IHC_annotation %>% 
  filter(`Cell type` == "cells in posterior")

IHC_annotation %>% 
  filter(`Cell type` == "sweat ducts")


# 118 celltyper
# 185 tissue + cell