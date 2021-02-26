
# Distance calculation function 

dist_calc <- function(df,comp,m,id=NULL) {
  start_time <- Sys.time()
  dist_matrix <-
    df  %>% as_tibble(rownames = "enssscg_id") %>%
    column_to_rownames("enssscg_id") %>%
    # head(10000) %>%
    select(1:comp) %>%
    get_dist(method = m)
  
  end_time <- Sys.time()
  total_time <- end_time - start_time
  
  if(!is.null(id)) {
    return (list(id = paste(id, m),
       distance = dist_matrix,
       time = total_time))
  }
  
  else {
    return (list(distance = dist_matrix,
                 time = total_time))
  }
}

# Clusteirng function

clust <- function(dist, k = 10, m, genes, id = NULL) {
  start_time <- Sys.time()
  
  if (m %in% c("kmeans", "pam", "clara", "fanny", "hclust", "agnes", "diana")) {
    res <- eclust(dist, FUNcluster = m, k = k, nboot = 500)
  }
  
  if (m == "fastkmeans") {
    centroids <- KMeans_arma(dist %>% as.matrix(), clusters = k,
                             n_iter = 10, seed_mode = "random_subset")
    res <- predict_KMeans(dist %>% as.matrix(), centroids)
    class(res) <- "numeric"
  }
  
  
  if (m == "louvain" | m == "leiden") {
    
    if (m == "louvain") {
      alg = 1
    }
    
    if (m == "leiden") {
      alg = 3
    }
        
    dist <- dist %>%
      as.matrix() %>%
      set_colnames(genes) %>%
      set_rownames(genes)%>%
      as.dist()
    
      louv <-
        CreateSeuratObject(assay = "Exp",
                           counts = t(data.frame(row.names = genes, count = c(1:22007))))
    neighbors <-
      FindNeighbors(
        dist,
        k.param = 20,
        return.neighbor = FALSE,
        compute.SNN = TRUE,
        prune.SNN = 1/15,
        nn.method = "annoy", #Distance metric for annoy. Options include: euclidean, cosine, manhattan, and hamming
        n.trees = 50,
        annoy.metric = "euclidean",####mirar
        nn.eps = 0,
        verbose = TRUE,
        force.recalc = FALSE,
        l2.norm = FALSE,
        cache.index = FALSE
      )
    
    louv@graphs$Exp_snn <- neighbors$snn
    
    louv <- FindClusters(louv, graph.name = "Exp_snn", resolution = 10, algorithm = alg)
    
    res <- as.numeric(as.character(louv@meta.data$Exp_snn_res.10))
    
  }
  
  if (m == "dbscan") {
    res_dbscan <- 
      dist %>%
      dbscan::dbscan(eps = 100, minPts = 5) 
    res <- as.numeric(res_dbscan$cluster)
    
  }
  
  if (m == "SOM") {
    res_som <- som(dist %>% as.matrix())
    res <- res_som$unit.classif
  }
  
  if (m == "fuzzyc") {
    res_fuzzy <- FKM(dist %>% as.matrix(), k = k)
    res <- as.numeric(res_fuzzy$clus[,1])
  }
  
  
  if (m == "fasthclust"){
    res <- flashClust(dist, method = "ward", members = NULL)
    res <- res %>% cutree(k) %>% dplyr::as_tibble(rownames = "gene") %>%
      column_to_rownames("gene")
    res$value = as.factor(res$value)
  }
  
  end_time <- Sys.time()
  
  total_time <- end_time - start_time
  print(total_time)

    # Add progress bar?
  
  res <- data.frame(gene = genes, value = res)
  
  if(!is.null(id)) {
    return(list(cluster = res,
                time = total_time,
                id = paste(id,m)))
  }
  else {
    return(list(cluster = res,
                time = total_time))
    }
  }

#load("data/processed/distances.Rdata")
#load("data/processed/gene_names.Rdata")

#d <- distances[[1]]
#rm(distances)

#graph <- graph.adjacency(d$distance %>% as.matrix(), mode = "undirected", weighted = TRUE, diag = TRUE)

#fromto <- get.edgelist(graph)
#clustergraph <- graph_from_edgelist(fromto, directed = FALSE) %>%
 # set.edge.attribute("weight", value = E(graph)$weight)

# Louvain parittion (weights)
#louvain_partition <- cluster_louvain(clustergraph, weights = E(graph)$weight)

# Community detection
#clustergraph$community <- louvain_partition$membership

#res3 <- clustergraph$community