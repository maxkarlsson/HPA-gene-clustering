
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
  # else {
  #   d <- dis[[1]]$distance
  # }
  
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
    scale_color_viridis(discrete = T) +
    theme(text = element_text(size=10)) + 
    labs(y = "Value", x = "Number of clusters")
  
  
  return(list(res = res,
              plot_k = plotk))
}



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

enrich_cluster <- function (clustering,num,ontology_type="BP",orthologs) {
  clustering2 <-
    clustering %>%
    filter(value == num) %>%
    left_join(orthologs, by = c("gene" = "enssscg_id")) %>%
    na.omit() %>%
    dplyr::select(entrez,value)
  
  gene_list <- clustering2$entrez
  # Profile clusters (with universe)
  ck <- try(enrichGO(gene = gene_list, 
                           OrgDb = org.Hs.eg.db, 
                           ont = ontology_type, 
                           universe = as.character(clustering$entrez),
                           pAdjustMethod = "BH",
                           pvalueCutoff = 1, 
                           qvalueCutoff = 1))
  
  ck_tibble <- as_tibble(ck) 

  
  return(ck_tibble)
  
}


find_bestk <- function (results) {
  #ank
  results <-
    results %>%
    mutate(rank_connectivity = rank(connectivity_index),
           rank_dunn = rank(-(dunn_index)),
           rank_asw = rank(-(avg_silhouette_width)),
           rank_BHI = rank(-(BHI_index)),
           rank_BPE = rank(-(BP_enriched)),
           rank_time = rank(time)
    )
  
  # Select best k 
  bestk <- 
    results %>%
    mutate(total_rank = rank_connectivity + rank_asw +  rank_BHI + rank_BPE +rank_dunn + rank_time) %>% #rank_BPE +rank_dunn + + rank_time
    dplyr::arrange(total_rank) %>%
    head(10)
  
    #group_by(scale_dist) %>%
    #mutate(min = (min(total_rank))) %>%
    #filter(total_rank == min)
  bestk$clustering_id <- factor(bestk$clustering_id,levels = bestk$clustering_id)
  
  plot_bestk <-
    bestk %>%
    as_tibble() %>%
    gather(key = measure, value = value, connectivity_index:time) %>%
    select(-scale_dist) %>%
    ggplot(aes(x=as.factor(clustering_id), y=value, fill=clustering_id)) +
    geom_bar(stat = "identity")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_wrap(~measure, scales="free_y") +
    scale_fill_viridis(discrete = T) +
    theme_minimal()+
    theme(axis.text.x=element_text(angle =- 50, vjust = 0.5)) + 
    scale_x_discrete(labels = NULL, breaks = NULL) +
    labs(y = "Value", x = "Clustering strategy")
  
  return(list(bestk= bestk,
              plot_bestk = plot_bestk))
  
}


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



