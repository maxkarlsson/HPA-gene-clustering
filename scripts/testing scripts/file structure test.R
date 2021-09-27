

library(tidyverse)



create_folder_structure <- 
  function(dataset_id, dataset_run_id, main_folder, folders) {
    dir.create(main_folder, showWarnings = FALSE)
    
    dataset_folder <-
      tibble(dataset_id,
             dataset_run_id,
             dataset_path = paste0(main_folder, "/",
                                   dataset_id, " ",
                                   dataset_run_id))
    
    for(path in dataset_folder$dataset_path) dir.create(path, showWarnings = FALSE)
    
    folders_df <-
      enframe(folders,
              "name", 
              "path") %>% 
      expand_grid(dataset_id) %>% 
      left_join(tibble(dataset_id,
                       dataset_run_id),
                by = "dataset_id") %>% 
      left_join(dataset_folder,
                by = c("dataset_id",
                       "dataset_run_id")) %>% 
      mutate(full_path = paste(dataset_path, path, sep = "/")) 
    
    for(path in folders_df$full_path) dir.create(path, showWarnings = FALSE)
    
    
    out_paths <- 
      set_names(dataset_folder$dataset_id,
                dataset_folder$dataset_id) %>% 
      lapply(function(x) {
        
        clustering_info <- 
          dataset_folder %>% 
          filter(dataset_id == x) 
        
        temp <- 
          folders_df %>% 
          filter(dataset_id == x) 
        
        list(main = clustering_info$dataset_path,
             run_id = clustering_info$dataset_run_id) %>% 
          append(split(temp$full_path,
                       temp$name))
        
      })
    
    
    return(out_paths)
  }

file_structure <- 
  create_folder_structure(dataset_id = c("tissue",
                                         "singlecell",
                                         "celline"),
                          dataset_run_id = c("test",
                                             "v2", 
                                             "v1"),
                          main_folder = "results/Clustering results",
                          folders = c("svg" = "svg",
                                      "treemap" = "svg/treemap",
                                      "heatmap" = "svg/heatmap",
                                      "UMAP" = "UMAP",
                                      "enrichment" = "enrichment",
                                      "clustering" = "clustering"))


file_structure[[dataset_id]]$clustering
