
# Create random clustering

random_cluster <- function (k,l, names) {
  start_time <- Sys.time()
  cluster <- sample(1:k, l, replace = T)
  end_time <- Sys.time()
  total_time <- end_time - start_time
  
  return(list(cluster = data.frame(gene = names, value = cluster),
              time = total_time,
              id = "random"))
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

eval <- function(clustering, genes, dis) {
  id <- clustering$id
  
  if(!grepl("2D",id)) {
    for (i in (1:length(dis))) {
      if(grepl(dis[[i]]$id,id)) {
        d <- dis[[i]]$distance
      }
    } 
  }
  else { #if (!grepl("2D",id))
    d <- dis$distance
  } 

  r <- 
    clustering$cluster %>%
    select(gene,value) %>%
    left_join(genes, by = c("gene" = "enssscg_id")) %>%
    na.omit() %>%
    dplyr::select(-gene) 
  
  clusters <- as.numeric(as.character(r$value))
  
  con_index <-
    connectivity(distance = d, cluster = as.numeric(as.character(clustering$cluster$value)))
  
  dunn_index <- 
    dunn(distance = d, cluster = as.numeric(as.character(clustering$cluster$value)))
  
  c <- silhouette(dist = d, 
                  x = as.numeric(as.character(clustering$cluster$value)))
  
  c_stats <- summary(c)
  
  asw_index <- c_stats$avg.width
  
  
  bio_index <- 
    if(require("Biobase") && require("annotate") && require("GO.db") &&
       require("org.Hs.eg.db")) {
      BHI(clusters, annotation="org.Hs.eg.db", names=r$entrez, category="BP") 
    }
  
  
  cluster_list <- c()
  num_clus <- tail(sort(r$value), n=1)
  num_genes <- length(r$value)
  for (i in c(0:num_clus)) {
    genes <- r %>%
      filter(value==i) %>%
      pull(entrez) %>%
      as.character()
    cluster_list[[as.character(i)]] <- genes
  }
  
  # Profile clusters (with universe)
  ck <- try(compareCluster(geneClusters = cluster_list, 
                           fun = "enrichGO", 
                           OrgDb = org.Hs.eg.db, 
                           ont = "BP", 
                           universe = as.character(r$entrez),
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.05))
  
  ck_tibble <- as_tibble(ck) 
  
  bio2_index <-
    length(unique(as.numeric(ck_tibble$Cluster)))/num_clus
  
  
  return(data.frame(clustering_id = id,
                    connectivity_index = con_index,
                    dunn_index = dunn_index,
                    avg_silhouette_width = asw_index,
                    BHI_index = bio_index,
                    BP_enriched = bio2_index,
                    time = as.numeric(clustering$time, units = "secs")
  ))
}


# Generate a data frame from evaluation reports
eval_to_df <- function (eval_list) {
  res_l <- c()
  res <- data.frame()
  
  for (i in c(1:length(eval_list))){
    res_l <- c(res_l,eval_list[[i]])
    res <- bind_rows(res,eval_list[[i]])
  }
  
  scaling <- c()
  distance <- c()
  clustering <- c()
  k <- c()
  
  for (row in 1:nrow(res)) {
    elements <- strsplit(res[[row, "clustering_id"]], " ")
    scaling <- c(scaling, elements[[1]][1])
    distance <- c(distance, elements[[1]][2])
    clustering <- c(clustering, elements[[1]][3])
    k <- c(k, elements[[1]][4])
  }
  
  res <-
    res %>%
    as_tibble %>%
    mutate(scaling = scaling,
           distance = distance,
           clustering = clustering,
           k = k, 
           scale_dist = paste(scaling,distance))
  
  to_plot <- res %>%
    as_tibble() %>%
    select(scale_dist,k,connectivity_index,dunn_index, avg_silhouette_width, 
           BHI_index, BP_enriched) %>% #time
    gather(key = measure, value = value, connectivity_index:BP_enriched) 
  
  plotk <- 
    to_plot %>%
    ggplot(aes(x=as.numeric(k), y=value, group =1, color = as.factor(measure))) +
    geom_point(show.legend = F) +
    geom_line(show.legend = F) +
    theme_minimal() +
    #ggtitle() +
    facet_wrap(scale_dist~measure, scales="free_y") +
    theme(text = element_text(size=10)) + 
    labs(y = "Value", x = "Number of clusters")
  
  
  return(list(res = res,
              plot_k = plotk))
}



# Calculate enrichment using clusterProfiler

enrich <- function (clustering,ontology_type="BP",orthologs) {
  clustering <-
    clustering %>%
    left_join(orthologs, by = c("gene" = "enssscg_id")) %>%
    na.omit() %>%
    dplyr::select(entrez,value)
  
  cluster_list <- c()
  num_clus <- tail(sort(clustering$value), n=1)
  num_genes <- length(clustering$value)
  for (i in c(0:num_clus)) {
    genes <- clustering %>%
      filter(value==i) %>%
      pull(entrez) %>%
      as.character()
    cluster_list[[as.character(i)]] <- genes
  }
  
  
# Calculate enrichments for a df
  
enrichments_df <- function(df, grouping, orthologs) {
    
    res <- df %>%
      group_by(.data[[grouping]]) %>%
      do({
        subset <- .
        
        enrichment <-
          enrich(clustering = subset, orthologs = orthologs) %>%
          mutate(id = unique(subset$id))
      })
    
  }
  
  
  # Profile clusters (with universe)
  ck <- try(compareCluster(geneClusters = cluster_list, 
                           fun = "enrichGO", 
                           OrgDb = org.Hs.eg.db, 
                           ont = ontology_type, 
                           universe = as.character(clustering$entrez),
                           pAdjustMethod = "BH",
                           pvalueCutoff = 1, 
                           qvalueCutoff = 1))
  
  ck_tibble <- as_tibble(ck) 
  
  bio_index <-
    length(unique(as.numeric(ck_tibble$Cluster)))/(num_clus+1)
  
  return(list(enriched_terms = ck_tibble,
              enrichment_score = bio_index))
  
}


# Calculate enrichment scores

enrichment_scores <- function(df,grouping) {
  
  df %>%
    group_by(.data[[grouping]]) %>%
    do ({
      data <- .
      min_significance <- data %>%
        mutate(log_qvalue = -log10(qvalue)) %>%
        group_by(Cluster) %>% 
        summarize(max(log_qvalue),na.rm = T) %>% 
        pull(`max(log_qvalue)`) %>% na.omit() %>% mean() 
      annotability <- length(data %>%
                               filter(qvalue < 0.05) %>% pull(Cluster) %>% sort() %>% unique())/
        length(data %>% pull(Cluster) %>% sort() %>% unique())
      data.frame(id = unique(data$id), min_significance = min_significance, annotability = annotability)
      
    })
}

enrichment_scores_kmeans <-
  map(enrichments_kmeans,
      ~enrichment_scores(df = .x, grouping = "id")) 



# Find clustering with the best k

find_bestk <- function (results,measures,plot = F, n=1) {
  best_id <-
    results %>%
    filter(measure %in% measures) %>%
    group_by(measure) %>%
    do ({
      data <- .
      if (unique(data$measure) == "connectivity_index") {
        data %>%
          mutate(rank= rank(value,ties.method = "average"))
      }
      else {
        data %>%
          mutate(rank = rank(-(value),ties.method = "average"))
      }
    }) %>%
    group_by(id) %>%
    summarise(total_rank = sum(rank)) %>%
    arrange(total_rank) %>% 
    head(n) %>%
    pull(id)
  
  if(plot == T) {
    results$k <- factor(results$k,levels = sort(unique(results$k)))
    
    title <- unique(results$scale_dist_clust)
    plot_bestk <-
      results %>%
      filter(measure %in% measures)  %>%
      ggplot(aes(x=k, y=value, fill=k)) +
      geom_bar(stat = "identity")+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~measure, scales="free_y") +
      scale_fill_viridis(discrete = T) +
      theme_minimal()+
      scale_x_discrete(labels = NULL, breaks = NULL) +
      labs(y = "Value", x = "Number of clusters") +
      ggtitle(title)
    
    
    return(list(bestk= best_id,
                plot_bestk = plot_bestk))
  }
  else {
    return(best_id)
  }
  
}

# Rank clustering 

clust_rank <- function (results,measures) {
  ranked <-
    results %>%
    filter(measure %in% measures) %>%
    group_by(measure) %>%
    do ({
      data <- .
      if (unique(data$measure) == "connectivity_index") {
        data %>%
          mutate(rank= rank(value,ties.method = "average"))
      }
      else {
        data %>%
          mutate(rank = rank(-(value),ties.method = "average"))
      }
    }) %>%
    group_by(id) %>%
    mutate(avg_rank = mean(rank)) %>%
    group_by(id) %>%
    arrange(avg_rank) 
  
  return(ranked)
}


# Add 2D and random clustering

add_clusterings <- function(df, ref_clust, dist, pca, gene_names, dist_m = "euclidean",
                            clust_m = "fasthclust", orthologs) {
  
  # 2D clustering
  gene_matrix <- pca$scores %>%
    set_rownames(gene_names) %>%
    as.data.frame() %>%
    select(c(1:30))
  
  umap <- gene_matrix %>%
    rownames_to_column("gene") %>%
    umap_calc(n_neigh = 15, row_names = "gene")
  
  dist_2d <- dist_calc(df = umap, m = dist_m, comp = 2)
  
  clustering_2d <- 
    clust(dist = dist_2d$distance, m = clust_m, r = 4, 
          genes = gene_names)
  
  clustering_2D <- list(cluster = clustering_2d$cluster, 
                        id = paste("zscore ", dist_m, "2D ", clust_m,  sep = ""), 
                        time = clustering_2d$time)
  
  eval_2D <- eval(clustering = clustering_2D, genes = ortholog_info, dis = dist_2d)
  
  #Random
  clustering_random <- ref_clust
  start_time <- Sys.time()
  clustering_random <-
    clustering_random %>%
    mutate(value = sample(value))
  end_time <- Sys.time()
  total_time <- end_time - start_time
  
  clustering_random = list(cluster = clustering_random, 
                           id = paste("zscore ", dist_m, "random ", clust_m,  sep = ""), 
                           time = total_time)
  
  eval_random <- eval(clustering = clustering_random, genes = ortholog_info, dis = dist)
  
  # Transform
  ev_extra <- list(eval_2D,eval_random)
  ev_extra <- eval_to_df(ev_extra)
  
  return(bind_rows(df,ev_extra$res))
  
}


# Transform lists of results into data frames

to_df <- function(multk, add_random = T) {
  
  lapply(c(1:4),
         function(i) 
           
           if(add_random == T) {
             random <- as.data.frame(multk[[1]][[i]]$cluster) %>%
               mutate(value = sample(value))
             df <- random %>%
               mutate(id = 
                        paste((gsub('[[:digit:]]+', '', 
                                    multk[[1]][[i]]$id)),"random", 
                              sep = ""))
             
             for (j in c(1:6)) {
               df <- bind_rows(df, as.data.frame(multk[[j]][[i]]$cluster) %>% 
                                 mutate(id = multk[[j]][[i]]$id))
             }
             return(df)
           }
         
         else {
           df <- data.frame()
           for (j in c(1:6)) {
             df <- bind_rows(df, as.data.frame(multk[[j]][[i]]$cluster) %>% 
                               mutate(id = multk[[j]][[i]]$id))
           }
           return(df)
         }
         
  )
}


# Merge data frames

merge_all_df <- function(list_df) {
  df <- data.frame()
  for (i in c(1:6)) {
    for (j in c(1:4)) {
      df <- bind_rows(df, list_df[[i]][[j]])
    }
  }
  
  return(df)
}




# Calculate omega squared from ANOVA results

omega_sq <- function(aov_in, neg2zero=T){
  aovtab <- summary(aov_in)[[1]]
  n_terms <- length(aovtab[["Sum Sq"]]) - 1
  output <- rep(-1, n_terms)
  SSr <- aovtab[["Sum Sq"]][n_terms + 1]
  MSr <- aovtab[["Mean Sq"]][n_terms + 1]
  SSt <- sum(aovtab[["Sum Sq"]])
  for(i in 1:n_terms){
    SSm <- aovtab[["Sum Sq"]][i]
    DFm <- aovtab[["Df"]][i]
    output[i] <- (SSm-DFm*MSr)/(SSt+MSr)
    if(neg2zero & output[i] < 0){output[i] <- 0}
  }
  output <- c(output, 1 - sum(output))
  names(output) <- c(rownames(aovtab)[1:n_terms], "Residuals")
  
  return(output)
} 



# Scale performance result and rank

scaled_rank <- function(results,measures) {
  results %>% 
    filter(measure %in% measures) %>%
    group_by(measure) %>%
    do ({
      data <- .
      if (unique(data$measure) == "connectivity_index") {
        data %>%
          mutate(scaled_value=  scales::rescale(-value,c(0,1)))
      }
      else {
        data %>%
          mutate(scaled_value=  scales::rescale(value,c(0,1)))
      }
    }) %>%
    ungroup() %>% 
    group_by(id) %>%
    mutate(avg_value = mean(scaled_value)) %>%
    group_by(id) %>%
    arrange(-avg_value)  
  
}


