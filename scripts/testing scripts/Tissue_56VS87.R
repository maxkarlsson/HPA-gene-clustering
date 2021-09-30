


library(tidyverse)
library(multidplyr)
library(clusterProfiler)
library(patchwork)
clustering <- 
  list(tissue56 = read_tsv("data/temp/tissue_56.tsv"),
     tissue87 = read_tsv("data/temp/tissue_87.tsv")) %>% 
  bind_rows(.id = "type")

tissue_class <- 
  read_tsv("../../Data/HPA/HPA21_E103_files/consensus_all_category_103.tsv")

tissue_class_terms <- 
  tissue_class %>% 
  mutate(enhanced_tissues = gsub(", ", ";", enhanced_tissues)) %>% 
  separate_rows(enhanced_tissues, sep = ",") %>% 
  mutate(enhanced_tissues = gsub(";", ", ", enhanced_tissues)) %>% 
  select(ensg_id, term = enhanced_tissues) %>% 
  filter(!is.na(term)) %>% 
  distinct()

  

# worker_cluster <- new_cluster(n = 2)
# cluster_library(worker_cluster, c("dplyr",
#                                   "tidyverse"))
# cluster_copy(worker_cluster, c("enricher"))
# 

term2gene <- 
  tissue_class_terms %>% 
  select(term, ensg_id)

# cluster_copy(worker_cluster, c("term2gene"))

# For 
res <- 
  clustering %>% 
  group_by(dataset_id, type) %>% 
  # partition(worker_cluster) %>% 
  do({
    g_clustering_data <<- .
    
    g_clustering_data %>% 
      group_by(cluster) %>% 
      
      do({
        g_data <- .
        
        pull(g_data, gene) %>% 
          enricher(maxGSSize = Inf, 
                   universe = unique(g_clustering_data$gene),
                   TERM2GENE = term2gene) %>% 
          as_tibble()
      }) %>% 
      ungroup()
    
  }) %>% 
  ungroup() %>% 
  collect()

# rm(worker_cluster)

plot_data_all <- 
  res %>% 
  mutate(GeneFrac = as.numeric(gsub("\\/.*", "", GeneRatio)) /
           as.numeric(gsub(".*\\/", "", GeneRatio)),
         BgFrac = as.numeric(gsub("\\/.*", "", BgRatio)) /
           as.numeric(gsub(".*\\/", "", BgRatio)),
         odds_ratio = GeneFrac / BgFrac) %>% 
  select(type, cluster, ID, odds_ratio)

plots <- 
  lapply(unique(plot_data_all$type), 
         function(id_) {
    
    plot_data <- 
      plot_data_all %>% 
      filter(type == id_)
    
    plot_clustering <- 
      plot_data %>% 
      select(cluster, ID, odds_ratio) %>% 
      spread(ID, odds_ratio, fill = 0) %>% 
      column_to_rownames("cluster") %>%
      list(tissue = t(.), 
           cluster = .) %>%
      map(. %>% 
            dist(method = "binary") %>% 
            hclust(method = "ward.D2"))
    
    plot_data %>%
      ungroup() %>% 
      mutate(cluster = factor(cluster,
                              with(plot_clustering$cluster,
                                   labels[order])),
             ID = factor(ID,
                         with(plot_clustering$tissue,
                              labels[order]))) %>%
      ggplot(aes(ID, cluster, 
                 size = odds_ratio, fill = odds_ratio)) +
      
      geom_point(shape = 21,
                 show.legend = T) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
      # coord_fixed() +
      scale_fill_gradient(low = "white", high = "orangered",
                          limits = c(0, max(plot_data_all$odds_ratio))) + 
      scale_size_continuous(range = c(0, 6), limits = c(0, max(plot_data_all$odds_ratio))) +
      ggtitle(id_)
  })


wrap_plots(plots)
ggsave("results/tissue 56 vs 87 bubbleheat.pdf",
       width = 12, height = 10)

plot_data_all %>% 
  ggplot(aes(type, odds_ratio)) +
  geom_violin(draw_quantiles = 0.5)
ggsave("results/tissue 56 vs 87 odds violin.pdf",
       width = 3, height = 5)


plot_data_all %>% 
  ggplot(aes(odds_ratio, fill = type)) +
  geom_density(alpha = 0.5)
ggsave("results/tissue 56 vs 87 odds density.pdf",
       width = 5, height = 3)

plot_data_all %>% 
  group_by(type) %>% 
  summarise(n = n_distinct(cluster))



