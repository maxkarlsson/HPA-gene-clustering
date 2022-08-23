

library(tidyverse)

infolder <- "results/Clustering_results/singlecell_HPA21v2/"
outfolder <- "results/Clustering_results/singlecell_HPA21v2_toIVIS/" 



if(dir.exists(infolder)) {
  files <- list.files(infolder, recursive = T)
}
  
    
if(!dir.exists(outfolder)) {
  dir.create(outfolder)
  
}
  
  
for(file in files) {
  
  
  folder <- 
    paste0(outfolder, file) %>% 
    str_extract(".*\\/")
  
  if(!dir.exists(folder)) {
    dir.create(folder, recursive = T)
  }
  
  
  
  if(str_detect(file, ".rds$")) {
    # browser()
    file_data <- read_rds(paste0(infolder, file))
    
    if(file == "graph/neighbors.rds") {
      new_file <- 
        paste0(outfolder, "graph/snn_neighbors.tsv")
      
    } else {
      new_file <- 
        str_replace(paste0(outfolder, file), "\\.rds$", ".tsv")
      
    }
    
    if(file.exists(new_file)) {
      cat(paste(file, "already exists. Moving on...\n"))
      next
    }
    
    # PCA
    if(file == "PCA/PCA.rds") {
      # browser()
      file_data$scores %>% 
        as_tibble(rownames = "gene") %>%
        write_tsv(new_file)
    }
    
    # Distance
    if(file == "distance/distances.rds") {
      # browser()
      file_data %>% 
        as.matrix() %>% 
        as_tibble(rownames = "gene") %>% 
        write_tsv(new_file)
    }    
    
    # Graph neighbors
    if(file == "graph/neighbors.rds") {
      # browser()
      file_data$snn %>%
        as.matrix() %>%
        as_tibble(rownames = "gene") %>%
        write_tsv(new_file)
    }
    
    
  } else {
    file.copy(paste0(infolder, file), 
              paste0(outfolder, file))
  }
    
    
    
}
