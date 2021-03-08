
# Create random clustering

random_cluster <- function (k,l, names) {
  cluster <- sample(1:k, l, replace = T)
  return(data.frame(gene = names, cluster = cluster))
}

# Transform results to dataframe
clusters_to_df <- function (res_list, gene_names) {
  df <- data.frame()
  times <- data.frame()
  
  for (cluster_method in res_list) {
    for (res in cluster_method) {
      new_entries <- data.frame(id = res$id, 
                                gene =  gene_names, 
                                cluster = as.numeric(as.character(res$cluster$value)))
      df <- bind_rows(df, new_entries)
      times <- bind_rows(times, data.frame(
                                id = res$id,
                                time = as.numeric(res$time, units = "secs")))
    }
  }
  return(list(cluster_results = df, cluster_time = times))
}

# Evaluation 

eval <- function(clustering, genes, dist, times) { # add distance # , t, id = NULL
  
  id <- unique(clustering$id)

  if (grepl("zscore euclidean",id)){d <- dist[[1]]$distance}
  if (grepl("min-max euclidean",id)){d <- dist[[2]]$distance}
  if (grepl("max euclidean",id)){d <- dist[[3]]$distance}
  if (grepl("zscore manhattan",id)){d <- dist[[4]]$distance}
  if (grepl("min-max manhattan",id)){d <- dist[[5]]$distance}
  if (grepl("max manhattan",id)){d <- dist[[6]]$distance}
  if (grepl("zscore pearson",id)){d <- dist[[7]]$distance}
  if (grepl("min-max pearson",id)){d <- dist[[8]]$distance}
  if (grepl("max pearson",id)){d <- dist[[9]]$distance}
  if (grepl("random",id)){d <- dist[[1]]$distance}
  
  
  r <- 
    clustering %>%
    select(gene,cluster) %>%
    left_join(genes, by = c("gene" = "enssscg_id")) %>%
    na.omit() %>%
    dplyr::select(-gene) 
  
  clusters <- as.numeric(as.character(r$cluster))
  
  con_index <-
    connectivity(distance = d, cluster = as.numeric(as.character(clustering$cluster)))
  
  dunn_index <- 
    dunn(distance = d, cluster = as.numeric(as.character(clustering$cluster)))
  
  bio_index <- 
    if(require("Biobase") && require("annotate") && require("GO.db") &&
       require("org.Hs.eg.db")) {
      BHI(clusters, annotation="org.Hs.eg.db", names=r$entrez, category="all") 
    }l
  
  t <- times %>%
    filter(id == id) %>%
    pull(time)
  
  #if(!is.null(unique(id))) {
    return(data.frame(clustering_id = id,
                      connectivity_index = con_index,
                      dunn_index = dunn_index,
                      BHI_index = bio_index,
                      time = t
                      ))
  #}
  
  #else {
  #  return(data.frame(connectivity = con_index,
  #                    dunn = dunn_index,
  #                    BHI = bio_index,
  #                    time = as.numeric(time, units = "secs")
  #                    ))
  #}

}
