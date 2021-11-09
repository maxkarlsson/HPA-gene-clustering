
# Connectivity & Average silhouette index


connectivity_file <- 
  paste0("data/processed/", run_id, "_connectivity.rds")
# silhouette_file <- 
#   paste0("data/processed/", run_id, "_silhouette.rds")

# Prepare cluster for parallelization
cluster <- new_cluster(n = 6)

cluster_library(cluster, c("dplyr",
                           "tidyverse",
                           "clValid",
                           "cluster"))

cluster_copy(cluster, c("connectivity_file",
                        "all_data_dists_filenames",
                        "all_data_cluster"))

Sys.time()

if(file.exists(connectivity_file)) {
  all_data_cluster_connectivity <- 
    readRDS(connectivity_file)
} else {
  
  
  all_data_cluster_connectivity <- 
    all_data_cluster %>% 
    group_by(dataset_id) %>% 
    partition(cluster) %>% 
    do({
      g_data <<- .
      g_dataset_id <- 
        unique(g_data$dataset_id)
      
      cat(paste("Running dataset:", g_dataset_id, "\n"))
      
      
      connectivity_dataset_file <-
        connectivity_file %>% 
        gsub(".rds$", "", .) %>% 
        paste0("_", g_dataset_id, ".rds") 
      
      
      if(file.exists(connectivity_dataset_file)) {
        cat(paste("  Files already present.\n"))
        
        g_res <- 
          readRDS(connectivity_dataset_file)
      } else {
        
        cat(paste("  Calculating metrics.\n"))
        
        
        temp_distan <- 
          readRDS(all_data_dists_filenames[g_dataset_id])
        
        g_res <- 
          g_data %>% 
          group_by(resolution, seed) %>% 
          summarise(connectivity = connectivity(distance = temp_distan, 
                                                cluster = cluster),
                    avg_silhouette_width = silhouette(dist = temp_distan, 
                                                      x = cluster) %>% 
                      summary() %>% 
                      {.$avg.width})
        # dunn = dunn(distance = temp_distan, 
        #             cluster = value))
        
        saveRDS(g_res, connectivity_dataset_file)
        
        
      }
      
      g_res
    }) %>% 
    collect() %>% 
    ungroup()
  
  
  saveRDS(all_data_cluster_connectivity, 
          connectivity_file)
}

Sys.time()


# if(file.exists(silhouette_file)) {
#   all_data_cluster_silhouette <- 
#     readRDS(silhouette_file)
# } else {
#   
#   all_data_cluster_silhouette <- 
#     all_data_cluster %>% 
#     group_by(dataset_id) %>% 
#     do({
#       g_data <<- .
#       
#       temp_distan <- 
#         readRDS(all_data_dists_filenames[unique(g_data$dataset_id)])
#       
#       
#       g_data %>% 
#         filter(resolution == 1, 
#                seed == 1) %>% 
#         group_by(resolution, seed) %>% 
#         summarise(silhouette = silhouette(dist = temp_distan, 
#                                           x = .$value) %>% 
#                     summary() %>% 
#                     {.$avg.width})
#           
#         })
#       
#       
#       a
#       
#       silhouette = silhouette(dist = d, 
#                               x = value)
#     })
#   
# 
#   
#   saveRDS(all_data_cluster_connectivity, 
#           silhouette_file)
# }


all_data_cluster_connectivity %>%
  gather(metric, value, connectivity, avg_silhouette_width) %>% 
  ggplot(aes(resolution, value, color = dataset_id)) +
  geom_point() +
  geom_smooth() +
  # scale_x_log10() +
  facet_grid(metric ~ dataset_id, scales = "free_y") +
  theme_bw() +
  ggtitle("Connectivity and Silhouette index")



all_data_cluster_connectivity %>% 
  left_join(all_data_cluster %>% 
              group_by(dataset_id,resolution, seed) %>% 
              mutate(k = max(cluster) + 1) %>% 
              select(-gene,-cluster) %>% 
              distinct) %>% 
  gather(metric, value, connectivity, avg_silhouette_width) %>% 
  ggplot(aes(k, value, color = as.factor(resolution))) +
  geom_point() +
  geom_smooth() +
  facet_grid(metric ~ dataset_id, scales = "free_y") +
  theme_bw() +
  ggtitle("Connectivity and Silhouette index")

# 
# ## Test - SAMPLE 
# MI_GOterms_randomized <- MI_GOterms %>%  mutate(ensg_id = sample(ensg_id))
# worker_cluster <- new_cluster(n = 6)
# cluster_library(worker_cluster, c("dplyr",
#                                   "tidyverse",
#                                   "infotheo"))
# cluster_copy(worker_cluster, c("MI_GOterms_randomized",
#                                "calculate_MI",
#                                "file_structure"))
# 
# 
# # Calculate mutual information (GO database)
# t <- Sys.time()
# 
# MI_res_consensus_randomized <-
#   all_clusters_consensus %>%
#   group_by(dataset_id) %>%
#   partition(worker_cluster) %>%
#   do({
#       g_data <- .
#     
#     g_dataset_id <- 
#       unique(g_data$dataset_id)
#     
#     MI_scores_randomized <- 
#         g_data %>% 
#         group_by(resolution) %>% 
#         
#         do({
#           g_data <- .
#           select(g_data, gene, cluster) %>%
#             calculate_MI(., MI_GOterms_randomized, nrep = 100)
#         }) %>% 
#         ungroup()
#       
#     MI_scores_randomized
# 
#   }) %>%
#   ungroup() %>%
#   collect()
# rm(worker_cluster)
# Sys.time() - t
# 
# 
# MI_res_consensus_randomized %>% 
#   left_join(res_k) %>% 
#   ggplot(aes(k,z_score, color = dataset_id)) +
#   geom_point() +
#   geom_smooth() +
#   theme_bw() +
#   scale_color_manual(values = dataset_palette) +
#   facet_grid(GO_domain~dataset_id, scales = "free_y")
# 
# MI_res_consensus %>% 
#   left_join(res_k) %>% 
#   ggplot(aes(k,z_score, color = dataset_id)) +
#   geom_point() +
#   geom_smooth() +
#   theme_bw() +
#   scale_color_manual(values = dataset_palette) +
#   facet_grid(GO_domain~dataset_id, scales = "free_y")


##### REACTOME

#9.25
# 
# g_data <- final_clustering
#  
# 
# MI_R <-
#   blood_selection %>% 
#   filter(seed == 1) %>% 
#   select(resolution, gene, cluster) %>%
#   mutate(cluster = as.character(cluster)) %>%
#   group_by(resolution) %>% 
#   do({
#       calculate_MI_reactome(., MI_Rterms, nrep = 300)
#   })
# 
# MI_GO <-
#   blood_selection %>% 
#   filter(seed == 1) %>% 
#   select(resolution, gene, cluster) %>%
#   mutate(cluster = as.character(cluster)) %>%
#   group_by(resolution) %>% 
#   do({
#       calculate_MI(., MI_GO, nrep = 300)
#   })
# 
# MI_R %>% 
#   ggplot(aes(resolution, z_score)) +
#   geom_point()
# 
# clust <- blood_selction_consenus %>% filter(resolution == "1.3")
# genes <- clust$gene
# 
# tissue_terms <- 
#   read_tsv("data/meta/consensus_all_category_92.tsv") %>% 
#   filter(ensg_id %in% genes) %>%
#   separate_rows("enhanced_tissues", sep = ',') %>%
#   select(ensg_id,enhanced_tissues) %>%
#   remove_missing() %>% 
#   rename(term_id = enhanced_tissues)
# 
# MI_tissue <-
#   blood_selection %>% 
#   filter(seed == 1) %>% 
#   select(resolution, gene, cluster) %>%
#   mutate(cluster = as.character(cluster)) %>%
#   group_by(resolution) %>% 
#   do({
#       calculate_MI_reactome(., tissue_terms, nrep = 300)
#   })
# 
#  calculate_MI_reactome(clust, tissue_terms, nrep = 300)
#  
#  MI_tissue %>% 
#    ggplot(aes(resolution, z_score)) +
#    geom_point()

