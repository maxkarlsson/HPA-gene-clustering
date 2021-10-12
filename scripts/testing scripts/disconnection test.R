
library(tidyverse)

graph_$snn %>% 
  dim

graph_$snn@i
graph_$snn@p
graph_$snn@simple_triplet_matrix

slam::simple_triplet_matrix

graph_$snn@p
graph_$snn@assay.used
graph_$snn@i
graph_$snn@Dim
graph_$snn@Dimnames
graph_$snn@x %>% length
graph_$snn@factors

graph_$snn[1, ] %>% 
  {.[which(. != 0)]} %>% 
  head(30)
graph_$snn@x %>% 
  head(30)

graph_tibble <- 
  graph_$snn %>% 
  as_tibble(rownames = "gene1")

graph_tibble_long <- 
  graph_tibble%>% 
  gather(gene2, edge, -1) %>% 
  filter(edge > 0)


A <- 
  graph_tibble_long %>% 
  group_by(gene1) %>% 
  mutate(n1 = n_distinct(gene2)) %>% 
  group_by(gene2) %>% 
  mutate(n2 = n_distinct(gene1)) %>% 
  group_by_all() %>% 
  mutate(n = min(c(n1, n2))) %>% 
  ungroup() %>% 
  arrange(n)

Apa <- 
  names(all_data_graphs) %>% 
  enframe("i",
          "dataset_id") %>% 
  group_by(dataset_id) %>% 
  do({
    g_data <- .
    
    graph_ <- all_data_graphs[[g_data$dataset_id]]
    
    graph_$snn %>% 
      as_tibble(rownames = "gene1") %>% 
      gather(gene2, edge, -1) %>% 
      filter(edge > 0) %>%  
      group_by(gene1) %>% 
      mutate(n1 = n_distinct(gene2)) %>% 
      group_by(gene2) %>% 
      mutate(n2 = n_distinct(gene1)) %>% 
      ungroup() 
  })

Apa

disconnected_genes <- 
  A %>% 
  filter(n == 1) %>% 
  pull(gene1)

all_clusters_consensus %>% 
  filter(gene %in% disconnected_genes,
         dataset_id == "singlecell",
         resolution == 4.7)

all_clusters_consensus %>% 
  filter(cluster == 17,
         dataset_id == "singlecell",
         resolution == 4.7)


all_data$singlecell %>% 
  filter(ensg_id %in% disconnected_genes) %>% 
  column_to_rownames("ensg_id") %>% 
  t() %>% 
  scale() %>% 
  pheatmap()

graph_$snn %>%
  {tibble(i = .@i, 
          p = .@p, 
          x = .@x)}

graph_$snn

########################################################
dim(graph_$snn)
UMAP_res <- 
  graph_$snn %>%
  # {.[1:nrow(.),1:nrow(.)]} %>% 
  RunUMAP(n.components = 2L, umap.method = "umap-learn") %>% 
  {.@cell.embeddings} %>% 
  as_tibble(rownames = "gene") 

graph_$snn %>%
  {.[1:nrow(.),1:nrow(.)]} %>% 
  class()

graph_$snn %>% 
  class()

UMAP_res %>% 
  left_join(all_clusters_consensus %>% 
              filter(resolution == 4.7,
                     dataset_id == "singlecell")) %>% 
  ggplot(aes(UMAP_1, UMAP_2, color = as.character(cluster))) +
  geom_point(size = 0.1,
             show.legend = F) +
  theme_void() +
  coord_fixed()




