
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
  if (m == "louvain") {
    res <- ""
  }
  if (m == "dbscan") {
    res <- ""
  }
  if (m == "SOM") {
    res <- som(dist %>% as.matrix())
  }
  if (m == "fuzzyc") {
    res <- ""
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

library(flashClust)
res <- flashClust(dist_to_cluster, method = "ward", members = NULL)
library(tidyverse)
res <- res %>% cutree(100) %>% as_tibble()
 