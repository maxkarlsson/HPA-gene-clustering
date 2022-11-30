library(tidyverse)
library(pheatmap)
library(ggplotify)
source("scripts/functions_utility.R")

dataset_metadata_21 <- 
  read_csv("run_settings/20221004_settings.csv")
            
dataset_metadata_22 <- 
  read_csv("run_settings/20220222_settings.csv")

data_paths_21 <- read_csv("run_settings/20211108 all_datasets.csv")
data_paths_22 <- read_csv("run_settings/20221020 all_datasets.csv")



file_structure_21 <-
  create_folder_structure(dataset_id = dataset_metadata_21$dataset_id,
                          dataset_run_id = dataset_metadata_21$version,
                          main_folder = "results/Clustering_results",
                          folders = c("clustering" = "clustering"))
file_structure_22 <-
  create_folder_structure(dataset_id = dataset_metadata_22$dataset_id,
                          dataset_run_id = dataset_metadata_22$version,
                          main_folder = "results/Clustering_results",
                          folders = c("clustering" = "clustering"))



clustering_data_21 <- 
  file_structure_21 %>% 
  map(function(x) x$clustering) %>% 
  map(function(x) {
    x <<- x
    file_ <- paste(x, "final_consensus.tsv", sep = "/")
    if(file.exists(file_)) {
      read_tsv(file_)
    } else {
      NULL
    }
  }) %>% 
  bind_rows(.id = "dataset_id") %>% 
  mutate(cluster = as.character(cluster)) 

clustering_data_22 <- 
  file_structure_22 %>% 
  map(function(x) x$clustering) %>% 
  map(function(x) {
    x <<- x
    file_ <- paste(x, "final_consensus.tsv", sep = "/")
    if(file.exists(file_)) {
      read_tsv(file_)
    } else {
      NULL
    }
  }) %>% 
  bind_rows(.id = "dataset_id") %>% 
  mutate(cluster = as.character(cluster)) 

####
list(v21 = clustering_data_21,
     v22 = clustering_data_22) %>% 
  bind_rows(.id = "version") %>% 
  group_by(dataset_id, version) %>% 
  count() %>% 
  ggplot(aes(n, dataset_id, fill = version, label = n)) +
  geom_col(position = position_dodge()) +
  geom_text(aes(group = version),
            position = position_dodge(width = 1), 
            hjust = 0) + 
  scale_x_continuous(expand = expansion(c(0, 0.5))) +
  theme_minimal()

###

joined_clusterings <- 
  clustering_data_21 %>% 
  left_join(clustering_data_22,
            by = c("dataset_id", "gene"),
            suffix = c("_21", "_22"))

overlap_res <- 
  joined_clusterings %>% 
  group_by(dataset_id, cluster_21) %>% 
  mutate(cluster_n_21 = n_distinct(gene)) %>% 
  group_by(dataset_id, cluster_22) %>% 
  mutate(cluster_n_22 = n_distinct(gene)) %>% 
  group_by(dataset_id, cluster_21, cluster_22, cluster_n_21, cluster_n_22) %>% 
  count() %>% 
  group_by_all() %>% 
  mutate(overlap = n / (cluster_n_21 + cluster_n_22 - n),
         overlap_max = n / max(c(cluster_n_21, cluster_n_22)))


overlap_res %>% 
  ungroup() %>% 
  filter(dataset_id %in% c("singlecell")) %>% 
  select(cluster_21,cluster_22,overlap_max) %>% 
  spread(cluster_22, overlap_max, fill = 0) %>% 
  column_to_rownames("cluster_21") %>% 
  pheatmap(clustering_method = "average",
           color = viridis::viridis(16),
           border_color = NA,
           cellheight = 9,
           cellwidth = 9
           ) %>% 
  as.ggplot()
ggsave(savepath("singlecell overlap heatmap singlecell.pdf"),
       width = 11, height = 11)
ggsave(savepath("singlecell overlap heatmap singlecell.png"),
       width = 11, height = 11)


overlap_res %>% 
  ungroup() %>% 
  filter(dataset_id == "singlecell") %>% 
  group_by(cluster_22) %>% 
  arrange(-overlap_max) %>% 
  slice(1) %>% 
  ungroup() %>% 
  # mutate(bin = cut(overlap_max, 16)) %>% 
  # mutate(bin_pos = bin %>% 
  #          str_remove("\\(") %>% 
  #          str_remove("\\]")) %>% 
  # 
  # separate(bin_pos, into = c("pos1", "pos2"), sep = ",") %>% 
  # mutate(pos1 = as.numeric(pos1),
  #        pos2 = as.numeric(pos2),
  #        center = (pos1 + pos2) / 2,
  #        width = (max(pos2) - min(pos1)) / 16) %>% %>% 
  
  ggplot(aes(overlap_max)) +
  geom_density(show.legend = F,
               fill = "gray90") +
  geom_tile(data = . %>% 
              summarise(min = min(overlap_max),
                        max = max(overlap_max),
                        width = (max - min) / (16),
                        bin = seq(min, max, length.out = 16) + width / 2),
            aes(x = bin, y = 0, fill = bin, height = width*4),
            show.legend = F) +
  scale_fill_viridis_c() +
  scale_y_continuous(expand = expansion(0)) +
  labs(x = "Overlap", y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
ggsave(savepath("cluster 21 v 22 max overlap distr.pdf"),
       width = 4, height = 3)
ggsave(savepath("cluster 21 v 22 max overlap distr.png"),
       width = 4, height = 3)

overlap_res %>% 
  ungroup() %>% 
  filter(dataset_id %in% c("celline")) %>% 
  select(cluster_21,cluster_22,overlap) %>% 
  spread(cluster_22, overlap, fill = 0) %>% 
  column_to_rownames("cluster_21") %>% 
  pheatmap(clustering_method = "ward.D2")


overlap_res %>% 
  ungroup() %>% 
  filter(dataset_id %in% c("singlecell")) %>% 
  select(cluster_21,cluster_22,overlap_max) %>% 
  spread(cluster_22, overlap_max, fill = 0) %>% 
  column_to_rownames("cluster_21") %>% 
  pheatmap(clustering_method = "ward.D2")

overlap_res %>% 
  ungroup() %>% 
  filter(dataset_id %in% c("celline")) %>% 
  select(cluster_21,cluster_22,overlap_max) %>% 
  spread(cluster_22, overlap_max, fill = 0) %>% 
  column_to_rownames("cluster_21") %>% 
  pheatmap(clustering_method = "ward.D2")




overlap_res %>% 
  filter(dataset_id %in% c("singlecell", "celline")) %>% 
  group_by(dataset_id, cluster_22) %>% 
  top_n(1, overlap) %>% 
  ggplot(aes(overlap, fill = dataset_id)) +
  geom_density(alpha = 0.5)


overlap_res %>% 
  
  group_by(dataset_id, cluster_22) %>% 
  mutate(max_overlap = max(overlap)) %>% 
  ungroup() %>% 
  filter(dataset_id %in% c("singlecell", "celline")) %>% 
  filter(overlap >= 0.1) %>% 
  arrange(desc(dataset_id), -max_overlap, -overlap) %>% 
  select(dataset_id, 
         old_cluster = 2,
         new_cluster = 3, 
         overlap) %>% 
  write_csv(savepath("v21 v22 cluster comparison.csv"))

