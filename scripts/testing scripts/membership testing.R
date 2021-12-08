
library(tidyverse)
library(readxl)

source("scripts/functions_utility.R")
source("scripts/cluster rasterizer.R")


dataset_metadata <- read_csv("run_settings/20210928_settings.csv")
file_structure <-
  create_folder_structure(dataset_id = dataset_metadata$dataset_id,
                          dataset_run_id = dataset_metadata$version,
                          main_folder = "results/Clustering_results",
                          folders = c("data" = "data", 
                                      "svg" = "svg",
                                      "UMAP" = "UMAP",
                                      "treemap" = "svg/treemap",
                                      "heatmap" = "svg/heatmap",
                                      "enrichment" = "enrichment",
                                      "clustering" = "clustering"))


membership_data <- 
  file_structure %>% 
  map(function(x) x$clustering) %>% 
  map(function(x) {
    file_ <- paste(x, "cluster_memberships.tsv", sep = "/")
    if(file.exists(file_)) {
      read_tsv(file_)
    } else {
      NULL
    }
  }) %>% 
  bind_rows(.id = "dataset_id") 

dataset_palette <- 
  c("tissue" = "#0083C4",  
    "brain" = "#FFDD00", 
    "blood" = "#CF161A", 
    "singlecell" = "#6AA692", 
    "celline" = "#97CF16")



membership_data %>% 
  filter(membership > 0.00001) %>% 
  mutate(membership = 1/membership) %>% 
  ggplot(aes(dataset_id, membership, fill = dataset_id)) +
  geom_violin(scale = "width",
              show.legend = F) +
  geom_text(data = . %>% 
              group_by(dataset_id) %>% 
              summarise(membership = max(membership)),
            aes(label = round(membership)), 
            vjust = -0.3) +
  scale_fill_manual(values = dataset_palette) +
  theme_void() +
  ggtitle("The Vases") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("doc/figures/The Vases.pdf",
       width = 4, height = 3)

membership_data %>% 
  group_by(gene) %>% 
  mutate(max_membership = max(membership)) %>% 
  filter(membership > 0.01) %>% 
  ggplot(aes(max_membership, membership, fill = dataset_id)) +
  geom_point(shape = 21, 
             show.legend = F) +
  scale_fill_manual(values = dataset_palette) +
  theme_void() +
  ggtitle("The Pyramid") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("doc/figures/The pyramid.pdf",
       width = 4, height = 3)

membership_data %>% 
  group_by(dataset_id, gene) %>% 
  mutate(max_membership = max(membership)) %>% 
  filter(membership > 0.01) %>% 
  ggplot(aes(max_membership, membership, fill = dataset_id)) +
  geom_point(shape = 21, 
             show.legend = F) +
  scale_fill_manual(values = dataset_palette) +
  theme_void() +
  ggtitle("The Flag") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("doc/figures/The Flag.pdf",
       width = 4, height = 3)

membership_data %>% 
  group_by(dataset_id, gene) %>% 
  summarise(max_membership = max(membership)) %>% 
  # filter(membership > 0.01) %>% 
  ggplot(aes(dataset_id, max_membership, fill = dataset_id)) +
  # geom_point(shape = 21, 
  #            show.legend = F) +
  geom_violin() +
  scale_fill_manual(values = dataset_palette) +
  # theme_void() +
  ggtitle("The Flag") +
  theme(plot.title = element_text(hjust = 0.5))

a <- 
  membership_data %>% 
  arrange(-membership) %>% 
  group_by(dataset_id, gene) %>% 
  # top_n(2, membership) %>% 
  slice(1:2) %>% 
  
  mutate(member_rank = paste("rank", row_number())) %>% 
  ungroup() 

View(a)

a %>% 
  select(-cluster) %>% 
  filter(gene == "ENSG00000000003") %>% 
  spread(member_rank, membership, fill = 0) %>% 
  # summarise(max_membership = max(membership),
  #           second_max_membership = membership[rank(membership, ties.method = "random") == 2]) %>% 
  # filter(max_membership < 0.5) %>% 
  # filter(membership > 0.01) %>% 
  ggplot(aes(`rank 1`, `rank 2`)) +
  # geom_point(shape = 21, 
  #            show.legend = F) +
  geom_hex(show.legend = F) +
  facet_wrap(~dataset_id) +
  # geom_violin() +
  # scale_fill_manual(values = dataset_palette) +
  theme_void() +
  ggtitle("The Rocks") +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_blank())
ggsave("doc/figures/The Rocks.pdf",
       width = 4, height = 3)


a %>% 
  select(-cluster) %>% 
  spread(member_rank, membership, fill = 0) %>% 
  # View
  # summarise(max_membership = max(membership),
  #           second_max_membership = membership[rank(membership, ties.method = "random") == 2]) %>% 
  # filter(max_membership < 0.5) %>% 
  # filter(membership > 0.01) %>% 
  ggplot(aes(`rank 1`, `rank 2`)) +
  # geom_point(shape = 21, 
  #            show.legend = F) +
  geom_hex(aes(fill = stat(log10(count))),
           show.legend = F) +
  facet_wrap(~dataset_id) +
  # geom_violin() +
  # scale_fill_manual(values = dataset_palette) +
  # theme_void() +
  ggtitle("The Rocks") +
  theme(plot.title = element_text(hjust = 0.5))


a %>% 
  select(-cluster) %>% 
  spread(member_rank, membership, fill = 0) %>% 
  # View
  # summarise(max_membership = max(membership),
  #           second_max_membership = membership[rank(membership, ties.method = "random") == 2]) %>% 
  # filter(max_membership < 0.5) %>% 
  # filter(membership > 0.01) %>% 
  ggplot(aes(`rank 1`, `rank 2`/`rank 1`)) +
  # geom_point(shape = 21, 
  #            show.legend = F) +
  geom_hex(aes(fill = stat(log10(count))),
           show.legend = F) +
  facet_wrap(~dataset_id) +
  # geom_violin() +
  # scale_fill_manual(values = dataset_palette) +
  # theme_void() +
  ggtitle("The Rocks") +
  theme(plot.title = element_text(hjust = 0.5))


membership_data %>% 
  
  group_by(dataset_id, gene) %>% 
  mutate(max_membership = max(membership)) %>% 
  arrange(max_membership) %>% 
  View


a %>% 
  select(-cluster) %>% 
  spread(member_rank, membership, fill = 0) %>% 
  mutate(show2nd = `rank 1` < 0.7 & `rank 2` >= 0.1) %>% 
  group_by(dataset_id, show2nd) %>% 
  count() 



#########################

cluster_annotation <-
  # list(tissue = read_csv("data/annotation/Tissue 20211013.csv"),
  #      blood = read_csv("data/annotation/Blood 20211013.csv"),
  #      celline = read_csv("data/annotation/Celline 20211013.csv"),
  #      brain = read_csv("data/annotation/Brain 20211013.csv")) %>% 
  list.files("data/annotation", full.names = T) %>% 
  {.[grepl("20211018", .)]} %>% 
  set_names(., gsub(".*\\d |\\.xlsx$", "", .)) %>% 
  map(. %>% read_excel()) %>% 
  bind_rows(.id = "dataset_id") %>% 
  rename(dataset_id = 1,
         cluster = 2,
         n_genes = 3,
         ann_specificity = 5, 
         ann_function = 6,
         reliability = 7,
         external_comment = 8,
         internal_comment = 9) %>% 
  mutate(cluster = as.character(cluster),
         full_annotation = paste(ann_specificity, 
                                 ann_function, 
                                 sep = " - "))


membership_data %>% 
  filter(dataset_id == "tissue") %>% 
  group_by(dataset_id, gene) %>% 
  filter(max(membership) < 0.6) %>% 
  ungroup() %>% 
  filter(membership > 0.1) %>%
  group_by(dataset_id, gene) %>% 
  filter(length(membership) > 1) %>% 
  ungroup() %>% 
  mutate(cluster = as.character(cluster)) %>% 
  left_join(cluster_annotation %>% 
              select(dataset_id, cluster, ann_specificity, ann_function)) %>% 
  arrange(dataset_id, gene, -membership) %>% 
  group_by(dataset_id, gene) %>% 
  slice(1, 2) %>% 
  View

