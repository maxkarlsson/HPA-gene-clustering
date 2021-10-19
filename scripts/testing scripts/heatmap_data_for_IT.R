





clustering_data %>% 
  filter(dataset_id == "singlecell") %>% 
  group_by(cluster) %>% 
  count() %>% 
  arrange(-n)
         cluster == 4)



heatmap_data <-
  readRDS("data/processed/combined_HPA_expression_data.RDS")


clustering_data %>% 
  filter(dataset_id == "singlecell") %>% 
  group_by(cluster) %>% 
  do({
    g_data <<- .
    
    plot_data <- 
      heatmap_data$singlecell_sample %>% 
      filter(ensg_id %in% g_data$gene) %>% 
      column_to_rownames("ensg_id") %>% 
      t() %>% 
      scale() %>% 
      t()
    
    
    plot_data[is.na(plot_data)] <- 0
    
    
    gene_hclust <- 
      plot_data %>% 
      dist() %>% 
      hclust(method = "ward.D2")
    
    sample_hclust <- 
      plot_data %>%
      t() %>% 
      dist() %>% 
      hclust(method = "ward.D2")
    
    plot_data[gene_hclust$order,
              sample_hclust$order] %>%
      as_tibble(rownames = "ensg_id") %>% 
      write_csv(paste0("results/Till Fredric/sc_cluster_",
                       unique(g_data$cluster),
                       ".csv"))
    
  })

clustering_data %>% 
  filter(dataset_id == "singlecell") %>% 
  group_by(cluster) %>% 
  do({
    g_data <<- .
    
    plot_data <- 
      heatmap_data$singlecell_consensus %>% 
      filter(ensg_id %in% g_data$gene) %>% 
      column_to_rownames("ensg_id") %>% 
      t() %>% 
      scale() %>% 
      t()
    
    
    plot_data[is.na(plot_data)] <- 0
    
    
    gene_hclust <- 
      plot_data %>% 
      dist() %>% 
      hclust(method = "ward.D2")
    
    sample_hclust <- 
      plot_data %>%
      t() %>% 
      dist() %>% 
      hclust(method = "ward.D2")
    
    plot_data[gene_hclust$order,
              sample_hclust$order] %>%
      as_tibble(rownames = "ensg_id") %>% 
      write_csv(paste0("results/Till Fredric/sc_consensus_",
                       unique(g_data$cluster),
                       ".csv"))
    
  })


heatmap_data$singlecell_sample %>% 
  filter()



worker_cluster <- new_cluster(n = 5)
cluster_library(worker_cluster, c("dplyr",
                                  "tidyverse"))
cluster_copy(worker_cluster, c("file_structure",
                               "cluster_membconf_data"))





heatmap_data_scaled %>% 
  group_by(dataset_id, heatmap_type, cluster) %>%  
  # partition(worker_cluster) %>%
  do({
    g_cluster_data <<- .
    
    g_dataset_id <- unique(g_cluster_data$dataset_id)
    
    g_heatmap_type <- unique(g_cluster_data$heatmap_type)
    
    
    plot_savename_zscore <- 
      paste0(file_structure[[g_dataset_id]]$heatmap, "/",
             g_heatmap_type,
             "_zscore_",
             unique(g_cluster_data$cluster), 
             ".svg")
    
    plot_savename_scaled <- 
      paste0(file_structure[[g_dataset_id]]$heatmap, "/",
             g_heatmap_type,
             "_scaled_",
             unique(g_cluster_data$cluster), 
             ".svg")
    
    
    g_cluster_data_wide <- 
      g_cluster_data %>% 
      select(ensg_id, sample, value_zscore) %>% 
      spread(sample, value_zscore) %>% 
      column_to_rownames("ensg_id") 
    
    gene_hclust <- 
      g_cluster_data_wide %>% 
      dist() %>% 
      hclust(method = "ward.D2")
    
    sample_hclust <- 
      g_cluster_data_wide %>%
      t() %>% 
      dist() %>% 
      hclust(method = "ward.D2")
    
    plot_data <- 
      g_cluster_data %>% 
      mutate(ensg_id = factor(ensg_id, 
                              rev(gene_hclust$labels[gene_hclust$order])),
             sample = factor(sample, 
                             sample_hclust$labels[sample_hclust$order])) 
    
    plot_confidence <- 
      cluster_membconf_data %>%
      filter(dataset_id == g_dataset_id) %>% 
      filter(ensg_id %in% gene_hclust$labels) %>% 
      filter(cluster == unique(g_cluster_data$cluster)) %>% #-> addition to filter memberships for other clusters
      mutate(ensg_id = factor(ensg_id, 
                              rev(gene_hclust$labels[gene_hclust$order]))) %>% 
      ggplot(aes("Confidence", ensg_id, fill = cluster_membership_confidence)) + 
      geom_tile_rast() +
      theme_void() +
      scale_fill_gradient(low = "#D1EEEA", high = "#2A5674",
                          limits = c(0.2, 1), # before limit to 0.5
                          breaks = c(0.5, 0.75, 1),
                          name = "Confidence") +
      theme(axis.text.x = element_text(angle = -90, face = "bold"), 
            legend.position = "left")
    
    plot_heat_zscore <- 
      plot_data %>% 
      ggplot(aes(sample, ensg_id, fill = value_zscore)) +
      geom_tile_rast() +
      scale_fill_gradient2(low = "#4575B4", high = "#D73027", mid = "#FFFFFF",
                           name = "Z-score") +
      theme_void() +
      theme(axis.text.x = element_text(angle = -90,
                                       hjust = 0, 
                                       size = 6),
            legend.position = "top")
    
    
    plot_heat_scaled <- 
      plot_data %>% 
      ggplot(aes(sample, ensg_id, fill = value_scaled)) +
      geom_tile_rast() +
      scale_fill_gradient2(low = "#4575B4", high = "#D73027", mid = "#FFFFFF",
                           name = "Relative Expression",
                           breaks = c(0, 0.5, 1),
                           limits = c(0, 1)) +
      theme_void() +
      theme(axis.text.x = element_text(angle = -90,
                                       hjust = 0,
                                       size = 6),
            legend.position = "top")
    
    column_width <- (1/53)
    plot_n <- length(sample_hclust$labels)
    plot_scaling_factor <- column_width * plot_n
    
    plot_widths <- 
      c(3, plot_n) %>% 
      {. * column_width} %>% 
      {. / max(.)}
    
    plot_width <- 
      1 + 4 * plot_scaling_factor
    
    p1 <- 
      plot_confidence + plot_heat_zscore +
      plot_layout(widths = plot_widths, nrow = 1)
    p2 <- 
      plot_confidence + plot_heat_scaled +
      plot_layout(widths = plot_widths, nrow = 1)
    
    ggsave(plot_savename_zscore, plot = p1,  
           width = plot_width, 
           height = 5)
    ggsave(plot_savename_scaled, plot = p2,
           width = plot_width, 
           height = 5)
    
    tibble()
  }) %>% 
  ungroup() %>% 
  collect()

rm(worker_cluster)