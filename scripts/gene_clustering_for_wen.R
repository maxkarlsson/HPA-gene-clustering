library(tidyverse)
library(magrittr)
library(uwot)
library(viridis)
library(pcaMethods)
library(factoextra)
library(Seurat)
library(uwot)


# ---- functions ----
scale_data <- function(df,
                       col_value,
                       col_gene,
                       col_sample, 
                       log_transform = F, 
                       scaling = "zscore") {
  
  if (log_transform == T) {
    df <-
      df %>%
      mutate(exp = log(.data[[col_value]] + 1))
  } else {
    df <-
      df %>%
      mutate(exp = .data[[col_value]])
  }
  
  if (scaling == "zscore") { 
    df_scaled <-
      df %>%
      group_by(.data[[col_gene]]) %>%
      mutate(mean = mean(exp),
             sd = sd(exp)) %>%
      ungroup() %>%
      mutate(exp = (exp - mean) / sd) 
  }
  
  else if (scaling == "min-max") {
    df_scaled <-
      df %>% 
      group_by(.data[[col_gene]]) %>%
      mutate(min = min(exp),
             max = max(exp)) %>%
      ungroup() %>%
      mutate (exp = (exp-min) / (max-min))
  }
  
  
  else if (scaling == "max") {
    df_scaled <-
      df %>% 
      group_by(.data[[col_gene]]) %>%
      mutate(max = max(exp)) %>%
      ungroup() %>%
      mutate (exp = exp/max) 
    
  }
  
  df_scaled %>% 
    select(.data[[col_gene]], .data[[col_sample]], exp) %>% 
    mutate(exp = ifelse(is.na(exp), 0, exp)) 
}

calculate_pca <- function(data, npcs) {
  require(pcaMethods)
  pca_res <-
    data %>%
    pca(nPcs = npcs)
  
  pca_stats <-
    summary(pca_res) %>%
    t() %>%
    as_tibble(rownames = "PC") %>%
    mutate(PC = as.numeric(gsub("PC", "", PC))) %>%
    dplyr::rename(R2cum = 3)
  
  list(pca = pca_res,
       scores = pcaMethods::scores(pca_res),
       loadings = pcaMethods::loadings(pca_res),
       stats = pca_stats) %>%
    return()
}

calculate_distance <- function(df, distance_metric) {
  
  df %>% 
    factoextra::get_dist(method = distance_metric)
  
}

cluster_genes <- function(genes,
                          neighbors,
                          clustering_method = "louvain",
                          k = 10, 
                          resolution = 1,
                          seed = seed) { 
  
  
  if (clustering_method == "louvain") {
    alg = 1
  } else if (clustering_method == "leiden") {
    alg = 4
  }
  
  louv <-
    CreateSeuratObject(assay = "Exp",
                       counts = t(data.frame(row.names = genes, c(1:length(genes)))))
  
  
  louv@graphs$Exp_snn <- neighbors$snn
  
  louv <- FindClusters(louv, graph.name = "Exp_snn", resolution = resolution, algorithm = alg,
                       random.seed = seed)
  
  col <- paste("Exp_snn_res.",resolution,sep="")
  res <- as.numeric(as.character(louv@meta.data[,col]))
  k <- tail(sort(as.numeric(as.character(louv@meta.data[,col]))), n=1)
  
  
  res <- data.frame(gene = genes, value = res)
  
  
  return(res)
}


generate_savefile <- 
  function(item, parameters, run_id, save_path = "data/processed/") {
    
    # Check if file exists
    paste0(save_path, run_id, "_", item, "_", parameters, ".rds")
    
  }



HPA_gene_clustering <- 
  function(data, 
           col_value,
           col_gene, 
           col_sample, 
           scaling = "zscore",
           distance_metric = "pearson",
           clustering_method = "louvain",
           npcs = 30,
           log_transform = F, 
           k = 10, 
           resolution = 1,
           seed = 42,
           run_id = "run") {
    
    # Scale data
    scaled_data <- 
      scale_data(df = data, 
                 col_value = col_value,
                 col_gene = col_gene, 
                 col_sample = col_sample,
                 log_transform = log_transform, 
                 scaling = scaling)
    
    # Perform dimensionality reduction using pca
    pca_data <- 
      scaled_data %>%
      spread(celltype, exp) %>% 
      column_to_rownames(col_gene) %>%
      calculate_pca(npcs = npcs)
    
    
    genes <- 
      rownames(pca_data$scores)
    
    # Calculate distance
    distance_data <- 
      pca_data$scores %>%
      calculate_distance(distance_metric = distance_metric) %>%
      as.matrix() %>%
      set_colnames(genes) %>%
      set_rownames(genes)%>%
      as.dist()
    
    # Calculate neighbor graph
    neighbors <-
      FindNeighbors(
        distance_data,
        k.param = 20,
        compute.SNN = TRUE,
        prune.SNN = 1/15,
        nn.method = "annoy", #Distance metric for annoy. Options include: euclidean, cosine, manhattan, and hamming
        annoy.metric = "euclidean",
        nn.eps = 0,
        verbose = TRUE,
        force.recalc = FALSE
      )
    
    
    
    # Cluster genes
    genes <- rownames(neighbors$nn)
    
    cluster_data <-
      crossing(resolution = resolution,
               seed = seed) %>% 
      group_by_all() %>% 
      do({
        cluster_genes(genes,
                      neighbors,
                      clustering_method = clustering_method,
                      resolution = .$resolution,
                      seed = .$seed)
      }) %>% 
      ungroup()
    
    
    
    
    list(pca_data = pca_data,
         neighbors = neighbors,
         cluster_data = cluster_data)
  }



#### calculate clusters -----


data <- 
  read_delim("data/consensus_hpa_92.tsv", delim = "\t") %>%
  select(ensg_id, celltype, tmm)

data_cluster <- 
  HPA_gene_clustering(data, 
                      col_value = "tmm",
                      col_gene = "ensg_id", 
                      col_sample = "celltype", 
                      scaling = "zscore",
                      distance_metric = "pearson",
                      clustering_method = "louvain",
                      seed = 1,
                      resolution = 1,
                      npcs = 30,
                      log_transform = F, 
                      run_id = "test")


cluster_palette <- 
  colorRampPalette(ggthemes::tableau_color_pal(palette = "Classic Cyclic")(13))(n_distinct(data_cluster$cluster_data$value))

  

umap_data <- 
  data_cluster$pca_data$scores %>%
  umap(n_neighbors = 15) %>% 
  as_tibble() %>% 
  mutate(ensg_id = rownames(data_cluster$pca_data$scores))

umap_data %>% 
  left_join(data_cluster$cluster_data,
            by = c("ensg_id" = "gene")) %>% 
  ggplot(aes(V1, V2, color = factor(value))) +
  geom_point(size = 0.1) +
  scale_color_manual(values = cluster_palette) +
  theme_bw()

