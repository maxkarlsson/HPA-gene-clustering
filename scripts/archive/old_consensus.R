# 
# ###########
# cluster <- new_cluster(n = 7)
# cluster_library(cluster, c("dplyr",
#                            "tidyverse",
#                            "magrittr",
#                            "clue")) 
# cluster_copy(cluster, c("find_consensus", 
#                         "to_cluster",
#                         "file_structure",
#                         "run_settings"))
# 
# # Find consensus clustering for each resolution
# t <- Sys.time()
# all_clusters_consensus <- 
#   all_data_cluster %>% 
#   group_by(dataset_id) %>% 
#   do({
#     g_data <- .
#     
#     g_dataset_id <- 
#       unique(g_data$dataset_id)
#     
#     savefile_cluster <- 
#       paste(file_structure[[g_dataset_id]]$clustering, "screening_consensus_2.tsv", sep = "/")
# 
#     
#     if(file.exists(savefile_cluster)) {
#       consensus <-
#         read_tsv(savefile_cluster)
#       
#     } else {
#       
#       consensus <- 
#         g_data %>% 
#         group_by(resolution) %>% 
#         partition(cluster) %>% 
#         do({
#           dat <- .
#           
#           final_clustering <- find_consensus(dat, n = run_settings$num_seeds, get_membership = F, runs = 1)
#           
#           final_clustering
#         }) %>% 
#         collect() %>% 
#         ungroup()
#       
#       write_tsv(consensus, savefile_cluster)
#       
#       consensus
#     }
#     
#   }) %>%  
#   ungroup()
# 
# rm(cluster)
# Sys.time() - t


```{r}
data_long <- 
  all_data_cluster %>%
  gather(seed, cluster, seed_1:seed_100) %>% 
  mutate(seed = as.numeric(gsub("seed_", "", seed)))

to_plot <- 
  data_long %>% 
  group_by(dataset_id, resolution, seed) %>% 
  summarise(k = n_distinct(cluster))

p1 <- to_plot %>% 
  select(dsataset_id, resolution, seed)
ggplot(aes(x = resolution, y = seed))  + 
  geom_bar(stat = "identity") +
  facet_grid(~dataset_id) +
  theme_bw() +
  scale_fill_manual(values = dataset_palette)

p2 <- data_long %>% 
  group_by(dataset_id, resolution, seed) %>% 
  summarise(k = n_distinct(cluster)) %>% 
  ggplot(aes(x = k, y = seed))  + 
  geom_bar(stat = "identity") +
  facet_grid(~dataset_id) +
  theme_bw() +
  scale_fill_manual(values = dataset_palette)

```



# cluster <- new_cluster(n = 5)
# cluster_library(cluster, c("dplyr",
#                            "tidyverse",
#                            "magrittr",
#                            "clue")) 
# cluster_copy(cluster, c("find_consensus", 
#                         "to_cluster",
#                         "run_settings",
#                         "file_structure"))
# 
# t <- Sys.time()
# # For all datasets
# final_clusterings_consensus <- 
#   final_clusterings %>% 
#   group_by(dataset_id) %>% 
#   partition(cluster) %>% 
#   do({
#     g_data <<- .
#     
#     g_dataset_id <- 
#       unique(g_data$dataset_id)
#     
#     savefile_cluster <- 
#       paste(file_structure[[g_dataset_id]]$clustering, "final_consensus.csv", sep = "/")
#     savefile_membership_matrix <- 
#       paste(file_structure[[g_dataset_id]]$clustering, "cluster_membership.csv", sep = "/")
#     savefile_cluster_mapping <- 
#       paste(file_structure[[g_dataset_id]]$clustering, "cluster_mapping_final_consensus.csv", sep = "/")
#     
#     if(file.exists(savefile_cluster)) {
#       consensus <-
#         read_csv(savefile_cluster)
#       
#     } else {
#       
#       consensus <- 
#         g_data %>% 
#         find_consensus(., id = "seed", runs = 5, get_membership = T) 
#       
#       write_csv(consensus$consensus_clustering, savefile_cluster)
#       write_csv(consensus$mapping_table, savefile_cluster_mapping)
#       write_csv(consensus$membership_matrix, savefile_membership_matrix)
#       
#       consensus <- consensus$consensus_clustering
#     }
#     consensus
#   }) %>% 
#   ungroup() %>% 
#   collect()
# 
# rm(cluster)
# Sys.time() - t
# 
# 
# final_clusterings_consensus %>% 
#   group_by(dataset_id,cluster) %>% 
#   summarize(a = n_distinct(gene)) %>% 
#   arrange(a)

#######################
# # Prepare cluster for parallelization
# cluster <- new_cluster(n = 6)
# 
# cluster_library(cluster, c("dplyr",
#                            "tidyverse",
#                            "factoextra",
#                            "pcaMethods",
#                            "Seurat",
#                            "magrittr",
#                            "cluster",
#                            "clValid")) 
# cluster_copy(cluster, c("HPA_gene_clustering", 
#                         "generate_savefile", 
#                         "scale_data",
#                         "calculate_pca",
#                         "calculate_distance",
#                         "cluster_genes", 
#                         "run_settings",
#                         "file_structure"))
# t <- Sys.time()
# 
# final_clusterings <- 
#   all_data %>% 
#   map(. %>% 
#         gather(sample_id, tmm, -1)) %>% 
#   bind_rows(.id = "dataset_id") %>% 
#   inner_join(expressed_genes) %>% 
#   group_by(dataset_id) %>%      
#   partition(cluster) %>% 
#   do({
#     g_data <- .
#     
#     g_dataset_id <- 
#       unique(g_data$dataset_id)
#     
#     savefile_cluster <- 
#       paste(file_structure[[g_dataset_id]]$clustering, "final_clustering.csv", sep = "/")
#     
#     if(file.exists(savefile_cluster)) {
#       cluster_res <-
#         read_csv(savefile_cluster)
#         
#     } else {
#       savefile_pca = paste(file_structure[[g_dataset_id]]$PCA, "PCA.rds", sep = "/")
#       savefile_pca_plot = paste(file_structure[[g_dataset_id]]$PCA, "PCA_plot.pdf", sep = "/")
#       savefile_dist = paste(file_structure[[g_dataset_id]]$distance, "distances.rds", sep = "/")
#       savefile_neighbors = paste(file_structure[[g_dataset_id]]$graph, "neighbors.rds", sep = "/")
#       
#       
#       cluster_res <-
#         HPA_gene_clustering(g_data,
#                             col_value = "tmm",
#                             col_gene = "ensg_id",
#                             col_sample = "sample_id",
#                             savefile_pca = savefile_pca,
#                             savefile_pca_plot = savefile_pca_plot,
#                             savefile_dist = savefile_dist,
#                             savefile_neighbors = savefile_neighbors,
#                             scaling = "zscore",
#                             distance_metric = "spearman",
#                             clustering_method = "louvain",
#                             seed = 1:100,
#                             resolution = run_settings$final_consensus_resolutions[[g_dataset_id]],
#                             npcs = 100,
#                             log_transform = F) 
#       
#       cluster_res <- 
#         cluster_res$cluster_data %>% 
#         mutate(cluster = as.character(cluster + 1))
#       
#       write_csv(cluster_res, savefile_cluster)
#     }
#     
#     cluster_res
#     
#   }) %>%
#   collect() %>% 
#   ungroup()
# 
# rm(cluster)
# Sys.time() - t