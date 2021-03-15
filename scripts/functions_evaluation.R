
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
  
  if ((id != "random")&(id != "2D")) {
    for (i in (1:length(dis))) {
      if(grepl(dis[[i]]$id,id)) {
        d <- dis[[i]]$distance
      }
    } 
  }
  else if (id == "2D"){
    d <- dis$distance
  } 
  else {
    d <- dis[[1]]$distance
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
  ck <- compareCluster(geneClusters = cluster_list, 
                             fun = "enrichGO", 
                             OrgDb = org.Hs.eg.db, 
                             ont = "BP", 
                             universe = as.character(r$entrez),
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05, 
                             qvalueCutoff = 0.05)
  
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

