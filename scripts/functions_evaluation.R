
# Create random clustering

random_cluster <- function (k,l, names) {
  cluster <- sample(1:k, l, replace = T)
  return(data.frame(gene = names, cluster = cluster))
}

# Transform results to dataframe
clusters_to_df <- function (res_list, gene_names) {
  df <- data.frame(id,gene,cluster)
  for res in res_list {
    new_entries <- data.frame(id = res$id, gene =  gene_names, cluster = res$cluster)
    df <- bind_rows(df, new_entries)
  }
  return(df)
}



# Evaluation 

eval <- function(clustering, genes, dist, id = NULL) { # add distance
  
  #clustering <- dplyr::rename(clustering, cluster = value)
  
  r <- 
    clustering %>%
    left_join(genes, by = c("gene" = "enssscg_id")) %>%
    na.omit() %>%
    dplyr::select(-gene) 
  
  
  clusters <- as.numeric(as.character(r$cluster))
  
  con_index <-
    connectivity(distance = dist, cluster = as.numeric(as.character(clustering$cluster)))
  
  dunn_index <- 
    dunn(distance = dist, cluster = as.numeric(as.character(clustering$cluster)))
  
  bio_index <- 
    if(require("Biobase") && require("annotate") && require("GO.db") &&
       require("org.Hs.eg.db")) {
      BHI(clusters, annotation="org.Hs.eg.db", names=r$entrez, category="all") 
    }

  
  if(!is.null(id)) {
    return(data.frame(connectivity = con_index,
                dunn = dunn_index,
                BHI = bio_index,
                id = id))
    
    #return(list(connectivity = con_index,
    #           dunn = dunn_index,
    #          BHI = bio_index,
    #         id = id))
  }
  else {
    return(data.frame(connectivity = con_index,
                dunn = dunn_index,
                BHI = bio_index))
  }

}