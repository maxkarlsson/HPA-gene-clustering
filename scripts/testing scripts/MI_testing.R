


library(tidyverse)
library(multidplyr)
library(infotheo)
library(clusterProfiler)
library(GOSemSim)

View(clusterProfiler::simplify)

clustering_data <-
  read_csv("data/processed/clusterings_to_annotate.csv") %>% 
  mutate(seed = 1)

GOterms <- 
  read_tsv("data/meta/Ensembl103 GO terms.txt") %>% 
  select(ensg_id = 1, 
         gene_name = 2, 
         GO_accession = 3, 
         GO_domain = 4, 
         GO_term_name = 5, 
         GO_term_evidence = 6) %>% 
  filter(ensg_id %in% clustering_data$gene)

clustering_data_seeds <- 
  readRDS("data/processed/clusters_5sds_5res.rds")

MI_GOterms <- 
  GOterms %>%
  filter(!is.na(GO_accession),
         !is.na(GO_domain)) %>%
  distinct() %>% 
  group_by(GO_accession) %>% 
  mutate(n_genes = n_distinct(ensg_id)) %>% 
  ungroup() %>% 
  filter(n_genes <= 500, 
         n_genes >= 3) %>% 
  select(ensg_id, GO_accession, GO_domain) %>% 
  distinct()



######

# sim <- mgoSim(res$ID, res$ID, semData = semData, measure = measure, 
#               combine = NULL)
# go1 <- go2 <- similarity <- NULL
# sim.df <- as.data.frame(sim)
# sim.df$go1 <- row.names(sim.df)
# sim.df <- gather(sim.df, go2, similarity, -go1)
# sim.df <- sim.df[!is.na(sim.df$similarity), ]
# sim.df <- merge(sim.df, res[, c("ID", by)], by.x = "go1", 
#                 by.y = "ID")
# sim.df$go2 <- as.character(sim.df$go2)
# ID <- res$ID
# GO_to_remove <- character()
# for (i in seq_along(ID)) {
#   ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > 
#                 cutoff)
#   if (length(ii) < 2) 
#     next
#   sim_subset <- sim.df[ii, ]
#   jj <- which(sim_subset[, by] == select_fun(sim_subset[, 
#                                                         by]))
#   GO_to_remove <- c(GO_to_remove, sim_subset$go1[-jj]) %>% 
#     unique
# }
# res[!res$ID %in% GO_to_remove, ]

###############################

# nrep = 100
# 
# 
# random_clustering <- 
#   clustering_data %>% 
#   group_by(dataset_id, resolution, seed) %>% 
#   do ({
#     dat <<- 
#       select(., gene, cluster)
#     
#     lapply(1:nrep,
#            function(x) {
#              dat %>% 
#                mutate(cluster = sample(cluster))
#            }) %>% 
#       bind_rows(.id = "run")
#            
#   }) %>% 
#   ungroup()


####

library(rbenchmark)

benchmark("original" = 
            {overlap_GO %>% 
                group_by(GO_domain) %>% 
                do({
                  domain_GOterms <- .
                  
                  if(random_simplify) {
                    domain_GOterms <- 
                      domain_GOterms %>% 
                      group_by(ensg_id) %>% 
                      summarise(GO_accession = sample(GO_accession, 1)) %>% 
                      ungroup()
                  }
                  
                  overlap_clustering %>%
                    expand_grid(run = 1:nrep) %>% 
                    group_by(run) %>% 
                    mutate(cluster = sample(cluster)) %>% 
                    ungroup() %>% 
                    inner_join(domain_GOterms,
                               by = c("gene" = "ensg_id")) %>% 
                    group_by(run) %>% 
                    summarise(MI = mutinformation(cluster, GO_accession))
                }) %>% 
                ungroup()},
          "separate_join" = 
            {overlap_GO %>% 
                group_by(GO_domain) %>% 
                do({
                  domain_GOterms <- .
                  
                  if(random_simplify) {
                    domain_GOterms <- 
                      domain_GOterms %>% 
                      group_by(ensg_id) %>% 
                      summarise(GO_accession = sample(GO_accession, 1)) %>% 
                      ungroup()
                  }
                  
                  expand_grid(run = 1:nrep) %>% 
                    group_by(run) %>% 
                    do({
                      overlap_clustering %>% 
                        mutate(cluster = sample(cluster)) %>% 
                        inner_join(domain_GOterms,
                                   by = c("gene" = "ensg_id")) %>% 
                        summarise(MI = mutinformation(cluster, GO_accession))
                    })
                    
                }) %>% 
                ungroup()},
          replications = 100)
          
          
t <- Sys.time()

Sys.time() - t

####

calculate_MI <-
  function(clustering, GO, nrep = 10, random_simplify = F) {
    overlap_GO <- 
      GO %>% 
      filter(ensg_id %in% clustering$gene)
    
    overlap_clustering <- 
      clustering %>% 
      filter(gene %in% overlap_GO$ensg_id)
    
    # t <- Sys.time()
    random_clustering_MI <- 
      overlap_GO %>% 
      group_by(GO_domain) %>% 
      do({
        domain_GOterms <- .
        
        if(random_simplify) {
          domain_GOterms <- 
            domain_GOterms %>% 
            group_by(ensg_id) %>% 
            summarise(GO_accession = sample(GO_accession, 1)) %>% 
            ungroup()
        }
        
        expand_grid(run = 1:nrep) %>% 
          group_by(run) %>% 
          do({
            overlap_clustering %>% 
              mutate(cluster = sample(cluster)) %>% 
              inner_join(domain_GOterms,
                         by = c("gene" = "ensg_id")) %>% 
              summarise(MI = mutinformation(cluster, GO_accession))
          })
      }) %>% 
      ungroup()
    # Sys.time() - t
    
    random_clustering_MI_specs <- 
      random_clustering_MI %>%
      group_by(GO_domain) %>% 
      summarise(mean_MI = mean(MI), 
                sd_MI = sd(MI))
    
    
    overlap_clustering_MI <- 
      overlap_GO %>% 
      group_by(GO_domain) %>% 
      do({
        domain_GOterms <- .
        
        if(random_simplify) {
          domain_GOterms <- 
            domain_GOterms %>% 
            group_by(ensg_id) %>% 
            summarise(GO_accession = sample(GO_accession, 1)) %>% 
            ungroup()
        }
        
        overlap_clustering %>%
          inner_join(domain_GOterms,
                     by = c("gene" = "ensg_id")) %>% 
          summarise(MI = mutinformation(cluster, GO_accession))
      }) %>% 
      ungroup() %>% 
      left_join(random_clustering_MI_specs,
                by = "GO_domain") %>% 
      mutate(z_score = (MI - mean_MI) / sd_MI)
    
    return(overlap_clustering_MI)
  }

nrep <- 20
MI_res <- 
  clustering_data %>% 
  # filter(dataset_id == "blood") %>% 
  group_by(dataset_id, resolution, seed) %>% 
  do({
    select(., gene, cluster) %>% 
      calculate_MI(., MI_GOterms, nrep = nrep)
  }) %>% 
  ungroup()

MI_res %>% 
  ggplot(aes(dataset_id, GO_domain, fill = z_score)) +
  geom_tile()

#####
t <- Sys.time()

worker_cluster <- new_cluster(n = 8)
cluster_library(worker_cluster, c("dplyr",
                                  "tidyverse",
                                  "infotheo"))
cluster_copy(worker_cluster, c("clustering_data_seeds",
                               "MI_GOterms",
                               "calculate_MI"))

MI_res_seeds <- 
  clustering_data_seeds %>% 
  group_by(dataset_id, resolution, seed) %>% 
  partition(worker_cluster) %>% 
  do({
    select(., gene, cluster) %>% 
      calculate_MI(., MI_GOterms, nrep = 50)
  }) %>% 
  ungroup() %>% 
  collect()

Sys.time() - t
#####

MI_res_seeds %>% 
  ggplot(aes(resolution, z_score)) +
  geom_line(data = . %>% 
              group_by(dataset_id, resolution, GO_domain) %>% 
              summarise(z_score = mean(z_score))) +
  geom_point() +
  facet_grid(GO_domain ~ dataset_id)

#####
t <- Sys.time()

worker_cluster <- new_cluster(n = 8)
cluster_library(worker_cluster, c("dplyr",
                                  "tidyverse",
                                  "infotheo"))
cluster_copy(worker_cluster, c("clustering_data_seeds",
                               "MI_GOterms",
                               "calculate_MI"))

MI_res_seeds_randsimp <- 
  clustering_data_seeds %>% 
  group_by(dataset_id, resolution, seed) %>% 
  partition(worker_cluster) %>% 
  do({
    select(., gene, cluster) %>% 
      calculate_MI(., MI_GOterms, nrep = 50, random_simplify = T)
  }) %>% 
  ungroup() %>% 
  collect()

Sys.time() - t
#####


MI_res_seeds_randsimp %>% 
  ggplot(aes(resolution, z_score)) +
  geom_line(data = . %>% 
              group_by(dataset_id, resolution, GO_domain) %>% 
              summarise(z_score = mean(z_score))) +
  geom_point() +
  facet_grid(GO_domain ~ dataset_id)

MI_res_seeds %>% 
  left_join(MI_res_seeds_randsimp,
            by = c("dataset_id", "resolution", "seed", "GO_domain"),
            suffix = c("", "randsimp")) %>% 
  ggplot(aes(z_score, z_scorerandsimp, color = resolution)) +
  geom_point() +
  facet_grid(GO_domain ~ dataset_id)

######################

t <- Sys.time()

worker_cluster <- new_cluster(n = 8)
cluster_library(worker_cluster, c("dplyr",
                                  "tidyverse",
                                  "infotheo"))
cluster_copy(worker_cluster, c("clustering_data_seeds",
                               "MI_GOterms",
                               "calculate_MI"))

MI_res_seeds_randsimp_iter1 <- 
  clustering_data_seeds %>% 
  filter(dataset_id == "tissue") %>% 
  expand_grid(run = c(1, 2),
              nrep = c(50, 100, 250)) %>% 
  group_by(dataset_id, resolution, seed, run, nrep) %>% 
  partition(worker_cluster) %>% 
  do({
    g_data <- .
    select(g_data, gene, cluster) %>% 
      calculate_MI(., MI_GOterms, nrep = unique(g_data$nrep), random_simplify = T)
  }) %>% 
  ungroup() %>% 
  collect()

MI_res_seeds_randsimp_iter2 <- 
  clustering_data_seeds %>% 
  filter(dataset_id == "tissue") %>% 
  expand_grid(run = c(1, 2),
              nrep = c(1000)) %>% 
  group_by(dataset_id, resolution, seed, run, nrep) %>% 
  partition(worker_cluster) %>% 
  do({
    g_data <- .
    select(g_data, gene, cluster) %>% 
      calculate_MI(., MI_GOterms, nrep = unique(g_data$nrep), random_simplify = T)
  }) %>% 
  ungroup() %>% 
  collect()

Sys.time() - t

MI_res_seeds_randsimp_iter <- 
  bind_rows(MI_res_seeds_randsimp_iter1, 
            MI_res_seeds_randsimp_iter2)

MI_res_seeds_randsimp_iter %>% 
  select(-contains("MI")) %>% 
  spread(run, z_score) %>% 
  ggplot(aes(`1`, `2`)) +
  geom_point() +
  geom_abline() +
  facet_grid(~nrep) +
  coord_fixed()



###################################################


t <- Sys.time()

worker_cluster <- new_cluster(n = 2)
cluster_library(worker_cluster, c("dplyr",
                                  "tidyverse",
                                  "infotheo"))
cluster_copy(worker_cluster, c("clustering_data_seeds",
                               "MI_GOterms",
                               "calculate_MI"))

MI_res_seeds_iter1 <- 
  clustering_data_seeds %>% 
  filter(dataset_id == "tissue") %>% 
  expand_grid(run = c(1, 2),
              nrep = c(50, 100, 250)) %>% 
  group_by(run, nrep, dataset_id, resolution, seed) %>% 
  partition(worker_cluster) %>% 
  do({
    g_data <- .
    select(g_data, gene, cluster) %>% 
      calculate_MI(., MI_GOterms, nrep = unique(g_data$nrep))
  }) %>% 
  ungroup() %>% 
  collect()

MI_res_seeds_iter2 <- 
  clustering_data_seeds %>% 
  filter(dataset_id == "tissue") %>% 
  expand_grid(run = c(1, 2),
              nrep = c(1000)) %>% 
  group_by(run, nrep, dataset_id, resolution, seed) %>% 
  partition(worker_cluster) %>% 
  do({
    g_data <- .
    select(g_data, gene, cluster) %>% 
      calculate_MI(., MI_GOterms, nrep = unique(g_data$nrep))
  }) %>% 
  ungroup() %>% 
  collect()

MI_res_seeds_iter <- 
  bind_rows(MI_res_seeds_iter1, 
            MI_res_seeds_iter2)

Sys.time() - t

MI_res_seeds_iter %>% 
  select(-contains("MI")) %>% 
  spread(run, z_score) %>% 
  ggplot(aes(`1`, `2`)) +
  geom_point() +
  geom_abline() +
  facet_grid(~nrep) +
  coord_fixed()


###################################################




t <- Sys.time()

worker_cluster <- new_cluster(n = 6)
cluster_library(worker_cluster, c("dplyr",
                                  "tidyverse",
                                  "infotheo"))
cluster_copy(worker_cluster, c("MI_GOterms",
                               "calculate_MI"))

MI_res_iter <- 
  clustering_data %>% 
  expand_grid(run = c(1, 2),
              nrep = c(10, 50, 100, 250, 1000)) %>% 
  group_by(run, nrep, dataset_id) %>% 
  partition(worker_cluster) %>% 
  do({
    g_data <- .
    select(g_data, gene, cluster) %>% 
      calculate_MI(., MI_GOterms, nrep = unique(g_data$nrep))
  }) %>% 
  ungroup() %>% 
  collect()

rm(worker_cluster)
gc()

Sys.time() - t

MI_res_iter %>% 
  select(-contains("MI")) %>% 
  spread(run, z_score) %>% 
  ggplot(aes(`1`, `2`)) +
  geom_point() +
  geom_abline() +
  facet_grid(~nrep) +
  coord_fixed()




#######


t <- Sys.time()

worker_cluster <- new_cluster(n = 6)
cluster_library(worker_cluster, c("dplyr",
                                  "tidyverse",
                                  "infotheo"))
cluster_copy(worker_cluster, c("MI_GOterms",
                               "calculate_MI"))

MI_res_iter <- 
  clustering_data %>% 
  group_by(dataset_id) %>% 
  partition(worker_cluster) %>% 
  do({
    g_data <- .
    select(g_data, gene, cluster) %>% 
      calculate_MI(., MI_GOterms, nrep = 1000)
  }) %>% 
  ungroup() %>% 
  collect()

rm(worker_cluster)
gc()

Sys.time() - t




##################################################

cluster <- new_cluster(n = 6)

cluster_library(cluster, c("dplyr",
                           "tidyverse",
                           "infotheo"))

cluster_copy(cluster, c("MI_file",
                        "all_data_cluster",
                        "rand_MI_data",
                        "GOterms"))

Sys.time()
if(file.exists(MI_file)) {
  all_data_cluster_MI <- 
    readRDS(MI_file)
} else {
  
  all_data_cluster_MI <-
    all_data_cluster %>% 
    group_by(dataset_id, resolution, seed) %>% 
    partition(cluster) %>% 
    do({
      g_data <<- .
      
      g_GOterms <- 
        GOterms %>% 
        filter(ensg_id %in% g_data$gene)
      
      g_data %>%
        ungroup() %>% 
        left_join(g_GOterms,
                  by = c("gene" = "ensg_id")) %>%
        summarise(MI = mutinformation(cluster, GO_accession))  
    }) %>% 
    collect() %>% 
    ungroup() %>%
    left_join(rand_MI_data) %>% 
    mutate(zscore = (MI - mean)/sd) 
  
  saveRDS(all_data_cluster_MI, 
          MI_file)
}
Sys.time()



all_data_cluster_MI %>%
  ggplot(aes(resolution, zscore, color = dataset_id)) +
  geom_point() +
  geom_smooth() +
  facet_grid( ~ dataset_id, scales = "free_y") +
  theme_bw() +
  ggtitle("Mutual information")

all_data_cluster_MI %>% 
  left_join(all_data_cluster %>% 
              group_by(dataset_id,resolution, seed) %>% 
              mutate(k = max(cluster) + 1) %>% 
              select(-gene,-cluster) %>% 
              distinct) %>% 
  ggplot(aes(resolution, zscore, color = as.factor(resolution))) +
  geom_point() +
  geom_smooth() +
  facet_grid(~dataset_id, scales = "free_y") +
  theme_bw() +
  ggtitle("Connectivity and Silhouette index")
