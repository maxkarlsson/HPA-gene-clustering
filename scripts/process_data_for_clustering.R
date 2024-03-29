

library(tidyverse)



data_path <- 
  "data/expression_data/HPA22_E103/"

dataset_path_file <- 
  "run_settings/20221020 all_datasets.csv"


dataset_paths <- 
  read_csv(dataset_path_file) %>% 
  mutate(full_path = paste0(data_path, 
                            filename))


##########

processing_functions <- 
  list(blood_sample = . %>% 
         select(ensg_id, sample, scilifelab_id, ntpm) %>% 
         filter(sample != "total PBMC") %>% 
         unite(sample, sample, scilifelab_id) %>% 
         spread(sample, ntpm),
       blood_consensus = . %>%
         select(ensg_id, sample = celltype, ntpm) %>%
         filter(sample != "total PBMC") %>% 
         spread(sample, ntpm),
       brain_tissue = . %>% 
         select(ensg_id, sample = tissue, ntpm) %>% 
         spread(sample, ntpm),
       brain_region = . %>% 
         select(ensg_id, sample = tissue, ntpm) %>%
         spread(sample, ntpm),
       tissue_region = . %>% 
         select(ensg_id, sample = tissue, ntpm) %>%
         spread(sample, ntpm),
       tissue_consensus = . %>% 
         select(ensg_id, sample = tissue, ntpm) %>%
         spread(sample, ntpm),
       celline_sample = . %>% 
         select(ensg_id, sample = cell_line, ntpm) %>%
         spread(sample, ntpm),
       celline_consensus = . %>% 
         select(ensg_id, sample = celline, ntpm) %>%
         spread(sample, ntpm),
       celline_disease = . %>% 
         select(ensg_id, sample = cancer, ntpm) %>%
         spread(sample, ntpm),
       singlecell_sample = . %>% 
         spread(sample, ntpm),
       singlecell_consensus = . %>% 
         spread(cell_type_name, exp))


##########
dataset_paths %>% 
  filter(!file.exists(full_path))

all_data <- 
  dataset_paths %>% 
  unite(id, dataset_id, type, sep = "_") %>% 
  filter(!is.na(filename)) %>%
  select(id, full_path) %>% 
  deframe() %>% 
  lapply(read_tsv)


processed_data <- 
  all_data %>%
  names() %>%
  set_names(., .) %>% 
  lapply(function(id) {
    
    one_dataset <- 
      all_data[[id]]
    
    dataset_processing_function <- 
      processing_functions[[id]]
    
    one_dataset %>% 
      dataset_processing_function()
  })


processed_data %>% 
  map(. %>% 
        select(-1) %>% 
        ncol())


saveRDS(processed_data, "data/processed/combined_HPA_expression_data.RDS")
saveRDS(processed_data, paste0("data/processed/", 
                               str_remove(dataset_path_file, ".*/") %>% 
                                 str_remove("\\.csv"), " combined_HPA_expression_data.RDS"))

