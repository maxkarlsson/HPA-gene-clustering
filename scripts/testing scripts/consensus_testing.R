
library(tidyverse)
library(vroom)
library(clue)
library(multidplyr)

screen_clusters <- 
  vroom("results/Clustering results/tissue testv2/clustering/screening_clustering.csv")


all_clusterings <- 
  screen_clusters# %>% 
  # filter(resolution %in% sample(unique(resolution), 5))


worker_cluster <- new_cluster(n = 5)
cluster_library(worker_cluster, c("dplyr",
                                  "tidyverse", 
                                  "clue"))
# cluster_copy(worker_cluster, c("cl_consensus"))

con_clusts <- 
  all_clusterings %>%
  filter(resolution > 5,
         resolution < 10) %>% 
  group_by(resolution) %>% 
  partition(worker_cluster) %>% 
  do({
    one_clustering <- .
    
    # Use the median cluster size as the consensus cluster size
    k <- 
      one_clustering %>% 
      select(id = seed, 
             cluster) %>% 
      group_by(id) %>% 
      summarise(n = n_distinct(cluster)) %>% 
      pull(n) %>% 
      median() %>% 
      ceiling()
    
    # Format clusters and calculate consensus (clue)
    clusters <-
      one_clustering %>% 
      select(id = seed, 
             gene,
             cluster) %>% 
      group_split(id) %>% 
      map(function(x) set_names(x$cluster, x$gene)) %>% 
      map(~as.cl_partition(.x))
    
    set.seed(1)
    cons_clustering <- cl_consensus(clusters, 
                                    method = "SE", 
                                    control = list(k = k, 
                                                   nruns = 100, 
                                                   verbose = FALSE)) 
    
    cons_clustering
    
    cl_class_ids(cons_clustering) %>% 
      enframe() %>% 
      mutate(cluster = as.character(value)) 
  }) %>% 
  group_by() %>% 
  collect()

rm(worker_cluster)
# final_clustering <- data.frame(gene = names(cl_class_ids(cons_clustering)), 
#                                cluster_cons = as.numeric(cl_class_ids(cons_clustering)),
#                                cluster = unclass(factor(as.numeric(cl_class_ids(cons_clustering)))))
con_clusts %>% 
  filter(!complete.cases(.))

con_clusts %>% 
  group_by(resolution, cluster) %>% 
  count() %>% 
  arrange(n, resolution) %>% 
  filter(n<5)

joined_seedcon_min <- 
  all_clusterings %>% 
  group_by(resolution, seed, cluster) %>% 
  count() %>% 
  group_by(resolution, seed) %>% 
  summarise(min_n = min(n)) %>% 
  left_join(con_clusts %>% 
              group_by(resolution, cluster) %>% 
              count() %>% 
              group_by(resolution) %>% 
              summarise(min_n = min(n)),
            by = "resolution", 
            suffix = c("_seed", "_con")) 

joined_seedcon_min %>% 
  arrange(min_n_con, min_n_seed) %>%
  filter(min_n_con < 5)
  View

joined_seedcon_min %>% 
  ggplot(aes(min_n_seed, min_n_con)) +
  geom_abline() +
  geom_point() +
  scale_y_continuous(limits = c(0, 50))
  

#######################################################

one_clustering <- 
  all_clusterings %>% 
  filter(resolution == 0.9)

# Use the median cluster size as the consensus cluster size
k <- 
  one_clustering %>% 
  select(id = seed, 
         cluster) %>% 
  group_by(id) %>% 
  summarise(n = n_distinct(cluster)) %>% 
  pull(n) %>% 
  median() %>% 
  ceiling()

# Format clusters and calculate consensus (clue)
clusters <-
  one_clustering %>% 
  select(id = seed, 
         gene,
         cluster) %>% 
  group_split(id) %>% 
  map(function(x) set_names(x$cluster, x$gene)) %>% 
  map(~as.cl_partition(.x))

set.seed(1)
cons_clustering <- cl_consensus(clusters, 
                                method = "SE", 
                                control = list(k = k, 
                                               nruns = 5, 
                                               verbose = FALSE)) 

cons_clustering

con_cluster <-
  cl_class_ids(cons_clustering) %>% 
  enframe() %>% 
  mutate(cluster = as.character(value)) 


con_cluster %>% 
  group_by(cluster) %>% 
  summarise(n_genes = n_distinct(name), 
            genes = paste(name, collapse = ";")) %>% 
  arrange(n_genes)

one_clustering %>% 
  filter(gene == "ENSG00000144161")

cons_clustering

cons_clustering$.Data["ENSG00000144161",] %>% 
  enframe %>% 
  ggplot(aes(as.character(name), value)) +
  geom_col()

cons_clustering$.Data["ENSG00000174891",] %>% 
  enframe %>% 
  ggplot(aes(as.character(name), value)) +
  geom_col()


gene_maxprob <- 
  cons_clustering$.Data %>% 
  apply(MARGIN = 1,
        function(x) max(x)) %>% 
  enframe("name", "max_prob")

cons_clustering$.Data[,12] %>% 
  enframe() %>% 
  left_join(con_cluster %>% 
              select(-value)) %>% 
  mutate(in_clust = cluster == 12) %>% 
  arrange(-value) %>% 
  left_join(gene_maxprob) %>% 
  mutate(ambig = value >= max_prob) %>% 
  filter(ambig)
