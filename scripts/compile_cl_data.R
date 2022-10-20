

library(tidyverse)


cl_data <- 
  read_tsv("../../Data/HPA/HPA22_E103_files/cell line v22/CCLE_2019_exp_103.tsv.gz")

cl_anno <- 
  readxl::read_excel("../../Data/HPA/HPA22_E103_files/cell line v22/CCLE_cell_line_annotation.xlsx")


cl_data_formatted <- 
  cl_data %>% 
  select(ensg_id, sample, ntpm) %>% 
  left_join(cl_anno %>% 
              select(sample = 1, 
                     cell_line = 4)) %>% 
  group_by(cell_line, ensg_id) %>% 
  summarise(ntpm = mean(ntpm))


sc_data_formatted %>% 
  write_tsv("data/expression_data/HPA22_E103/celline_CCLE_data.tsv")


