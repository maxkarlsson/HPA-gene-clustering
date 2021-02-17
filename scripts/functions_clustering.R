library(tidyverse)
library(ClusterR)
library(kohonen)
library(fclust)
library(flashClust)

# Distance calculation function 

dist_calc <- function(df,comp,m) {
  dist_matrix <-
    df  %>% as_tibble(rownames = "enssscg_id") %>%
    column_to_rownames("enssscg_id") %>%
   # head(10000) %>%
    select(1:comp) %>%
    get_dist(method = m)
}


# Clusteirng function

clust <- function(dist, k = 10, m) {
  start_time <- Sys.time()
  
  if (m %in% c("kmeans", "pam", "clara", "fanny", "hclust", "agnes", "diana")) {
    res <- eclust(dist, FUNcluster = m, k = k, nboot = 500)
  }
  
  if (m == "fastkmeans")) {
    centroids <- KMeans_arma(dist %>% as.matrix(), clusters = k,
                             n_iter = 10, seed_mode = "random_subset")
    res <- predict_KMeans(dist %>% as.matrix(), centroids, threads = 2)
  }
  
  if (m == "louvain") {
    # not tried - try with 500 genes first
    
    # Create graph from distance matrix
    graph <- grap.adjacency(dist %>% as.matrix(), mode = "undirected", weighted = TRUE. diag = TRUE)
    fromto <- get.edgelist(graph)
    clustergraph <- graph_from_edgelist(fromto, directed = FALSE) %>%
      set.edge.attribute("weight", value = E(graph)$weight)
    
    # Louvain parittion (weights)
    louvain_partition <- cluster_louvain(clustergraph, weights = E(graph)$weight)
    
    # Community detection
    clustergraph$community <- louvain_partition$membership
    
    res <- clustergraph$community
  }
  if (m == "dbscan") {
    #test
    #ADD library(dbscan) to main
    
    res <- 
      dist %>%
      dbscan::dbscan(eps = 10) #???
  }
  if (m == "SOM") {
    res <- som(dist %>% as.matrix())
  }
  if (m == "fuzzyc") {
    res <- FKM(dist %>% as.matrix(), k = k)
  }
  
  if (m == "fasthclust"){
    res <- flashClust(dist, method = "ward", members = NULL)
    res <- res %>% cutree(k) %>% dplyr::as_tibble(rownames = "gene") %>%
      column_to_rownames("gene")
    res$value = as.factor(res$value)
    
    pdf("dendogram.pdf")
    dist %>%
      pheatmap(annotation_row = res)

    res %>%
      fviz_dend(cex = 0.5, k = k)
    dev.off()
  }
  end_time <- Sys.time()
  
  total_time <- end_time - start_time
  print(total_time)
  # Add progress bar?
  return(res)
}

load("data/processed/distance_zscore.Rdata")


res <- flashClust(dist_to_cluster, method = "ward", members = NULL)
res <- res %>% cutree(100) %>% as_tibble()
 