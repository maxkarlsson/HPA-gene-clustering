
# ---- functions ----

elongate_sparsemat <- 
  function(x) {
    x_ <- 
      x %>% 
      summary() %>% 
      as_tibble()
    
    id1 <- rownames(x)[x_$i]
    id2 <- colnames(x)[x_$j]
    
    x_ %>% 
      mutate(id1 = id1,
             id2 = id2) %>% 
      select(id1, id2, x)
  }

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
       sdev = pcaMethods::sDev(pca_res),
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
  
  suppressWarnings(
    {
      louv <-
        CreateSeuratObject(assay = "Exp",
                           counts = t(data.frame(row.names = genes, c(1:length(genes)))), )
      louv@graphs$Exp_snn <- neighbors$snn
      
      louv <- FindClusters(louv, graph.name = "Exp_snn", resolution = resolution, algorithm = alg,
                           random.seed = seed, verbose = F)
    }
  )
  
  
  
  
  col <- paste("Exp_snn_res.",resolution,sep="")
  res <- as.numeric(as.character(louv@meta.data[,col]))
  k <- tail(sort(as.numeric(as.character(louv@meta.data[,col]))), n=1)
  
  
  res <- data.frame(gene = genes, cluster = res)
  
  
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
           savefile_pca,
           savefile_pca_plot,
           savefile_dist,
           savefile_neighbors,
           scaling = "zscore",
           distance_metric = "pearson",
           clustering_method = "louvain",
           npcs = 30,
           log_transform = F, 
           k = 10, 
           resolution = 1,
           seed = 42) {
    
    
    if(file.exists(savefile_pca)) {
      pca_data <- read_rds(savefile_pca)  
    } else {
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
        spread(col_sample, exp) %>% 
        column_to_rownames(col_gene) %>%
        calculate_pca(npcs = npcs)
      
      
      saveRDS(pca_data, savefile_pca)
    }
    
    ncomp <- # Kaiser's rule: retain factors whose eigenvalues are greater than 1 (they explain at least as much variance as one sample in the original dataset)
      pca_data$sdev[(pca_data$sdev)^2 >= 1] %>% # under some assumptions, variance = eigenvalues
      tail(1) %>%  
      enframe () %$% 
      name 
    ncomp <- as.numeric(gsub("PC", "", ncomp))
    
    if (pca_data$stats[ncomp, ]$R2cum < 0.8) {
      ncomp <- pca_data$stats %>% 
        filter(R2cum > 0.8) %>% 
        pull(PC) %>% 
        head(1)
    }

    plot <- 
      pca_data$stats %>%
      ggplot(aes(PC,R2cum)) +
      geom_point() +
      geom_line() +
      theme_bw() +
      theme_bw() +
      geom_vline(xintercept = ncomp, linetype = "dashed") +
      annotate("text",
               x = ncomp,
               y = 0.55,
               label = paste0("PC ", ncomp,
                              "\nR2 = ", round(pca_data$stats[ncomp, ]$R2cum, 3)),
               hjust = 1,
               vjust = 0)
    
    ggsave(savefile_pca_plot, plot = plot, width = 5, height = 5)
    
    
    if(file.exists(savefile_dist)) {
      distance_data <- read_rds(savefile_dist)
    } else {
      
      
      genes <- rownames(pca_data$scores)
      # Calculate distance
      distance_data <- 
        pca_data$scores %>%
        as.data.frame() %>% 
        select(1:ncomp) %>% # Select the number of components determined with the Kaiser rule
        calculate_distance(distance_metric = distance_metric) %>%
        as.matrix() %>%
        set_colnames(genes) %>%
        set_rownames(genes)%>%
        as.dist()
      
      saveRDS(distance_data, savefile_dist)
      
    }
    
    
    if(file.exists(savefile_neighbors)) {
      neighbors <- read_rds(savefile_neighbors)
    } else {
      
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
      
      saveRDS(neighbors, savefile_neighbors)
    }
    
    
    # Cluster genes
    genes <- rownames(neighbors$nn)
    
    cluster_data <-
      tidyr::crossing(resolution = resolution,
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


HPA_gene_clustering_old <- 
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
    
    ####
    savefile_pca <- 
      generate_savefile("pca", 
                        parameters = paste(scaling,
                                           npcs,
                                           log_transform,
                                           sep = "_"), 
                        run_id = run_id)
    
    savefile_dist <- 
      generate_savefile("dist", 
                        parameters = paste(scaling,
                                           distance_metric,
                                           npcs,
                                           log_transform,
                                           sep = "_"), 
                        run_id = run_id)
    
    savefile_neighbors <- 
      generate_savefile("neighbors", 
                        parameters = paste(scaling,
                                           distance_metric,
                                           npcs,
                                           log_transform,
                                           sep = "_"), 
                        run_id = run_id)
    ####
    
    if(file.exists(savefile_pca)) {
      pca_data <- read_rds(savefile_pca)  
    } else {
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
        spread(col_sample, exp) %>% 
        column_to_rownames(col_gene) %>%
        calculate_pca(npcs = npcs)
      
      
      saveRDS(pca_data, savefile_pca)
    }
    
    ncomp <- # Kaiser's rule: retain factors whose eigenvalues are greater than 1 (they explain at least as much variance as one sample in the original dataset)
      pca_data$sdev[(pca_data$sdev)^2 >= 1] %>% # under some assumptions, variance = eigenvalues
      tail(1) %>%  
      enframe () %$% 
      name 
    ncomp <- as.numeric(gsub("PC", "", ncomp))
    
    if (pca_data$stats[ncomp, ]$R2cum < 0.8) {
      ncomp <- pca_data$stats %>% 
        filter(R2cum > 0.8) %>% 
        pull(PC) %>% 
        head(1)
    }
    
    pca_data$stats %>%
      ggplot(aes(PC,R2cum)) +
      geom_point() +
      geom_line() +
      theme_bw() +
      theme_bw() +
      geom_vline(xintercept = ncomp, linetype = "dashed") +
      annotate("text",
               x = ncomp,
               y = 0.55,
               label = paste0("PC ", ncomp,
                              "\nR2 = ", round(pca_data$stats[ncomp, ]$R2cum, 3)),
               hjust = 1,
               vjust = 0) +
      ggtitle(run_id)
    
    ggsave(paste0("results/PCA_plot_",run_id,".pdf"))
    
    
    if(file.exists(savefile_dist)) {
      distance_data <- read_rds(savefile_dist)
    } else {
      
      
      genes <- rownames(pca_data$scores)
      # Calculate distance
      distance_data <- 
        pca_data$scores %>%
        as.data.frame() %>% 
        select(1:ncomp) %>% # Select the number of components determined with the Kaiser rule
        calculate_distance(distance_metric = distance_metric) %>%
        as.matrix() %>%
        set_colnames(genes) %>%
        set_rownames(genes)%>%
        as.dist()
      
      saveRDS(distance_data, savefile_dist)
      
    }
    
    
    if(file.exists(savefile_neighbors)) {
      neighbors <- read_rds(savefile_neighbors)
    } else {
      
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
      
      saveRDS(neighbors, savefile_neighbors)
    }
    
    
    # Cluster genes
    genes <- rownames(neighbors$nn)
    
    cluster_data <-
      tidyr::crossing(resolution = resolution,
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

to_cluster <- function(cluster_df) {
  v <- cluster_df$cluster
  names(v) <- cluster_df$gene
  return(v)
}

df_to_vector <- function(cluster_df) {
  v <- cluster_df$cluster
  names(v) <- cluster_df$gene
  return(v)
}


find_consensus <- function(all_clusterings, id, get_membership = F, runs = 1) {
  
  # Use the median cluster size as the consensus cluster size
  k <- 
    all_clusterings %>% 
    select(id = id, 
           cluster) %>% 
    group_by(id) %>% 
    summarise(n = n_distinct(cluster)) %>% 
    pull(n) %>% 
    median() %>% 
    ceiling()
  
  # Format clusters and calculate consensus (clue)
  clusters <- 
    all_clusterings %>% 
    select(id = id, 
           gene,
           cluster) %>% 
    group_split(id) %>% 
    map(to_cluster) %>%  
    map(~as.cl_partition(.x))
  
  set.seed(1)
  cons_clustering <- cl_consensus(clusters, 
                                  method = "SE", 
                                  control = list(k = k, 
                                                 nruns = runs, 
                                                 verbose = FALSE)) 
  
  cons_clusters <- as.numeric(cl_class_ids(cons_clustering)) %>% unique() %>% sort()
  final_clusters <- c(1:length(cons_clusters))
                            
             
  mapping_table <-  data.frame(cluster_cons = cons_clusters,
                               cluster = final_clusters)
  
  final_clustering <- data.frame(gene = names(cl_class_ids(cons_clustering)), 
                                 cluster_cons = as.numeric(cl_class_ids(cons_clustering))) %>% 
    left_join(mapping_table) %>% 
    select(-cluster_cons)
  
  
  if (get_membership) {
    # Extract cluster membership matrix 
    cons_matrix <- 
      cons_clustering$.Data[,] %>% 
      as_tibble(rownames = "gene") %>% 
      gather(cluster, membership, -1) %>% 
      filter(membership > 0) %>% 
      mutate(cluster = as.numeric(gsub("V", "", cluster))) %>% 
      rename(cluster_cons = cluster) %>% 
      left_join(mapping_table) %>% 
      select(-cluster_cons)
    
    return(list(consensus_clustering = final_clustering, 
                mapping_table = mapping_table,
                membership_matrix = cons_matrix ))
  } else {
    return(list(consensus_clustering = final_clustering, 
                mapping_table = mapping_table))
  }
}
