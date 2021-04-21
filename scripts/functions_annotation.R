plot_treemap <- 
  function(GO_terms,GO_ontology, clus) {
    plot_data <-
      GO_terms$enriched_terms %>% 
      filter(Cluster == clus) %>%
      filter(qvalue < 0.01)
    
    if(nrow(plot_data) > 1) {
      plot_mat <-
        calculateSimMatrix(plot_data$ID,
                           orgdb="org.Hs.eg.db",
                           ont=GO_ontology,
                           method="Rel")
      plot_mat2 <-
        reduceSimMatrix(plot_mat,
                        setNames(-log10(plot_data$qvalue), plot_data$ID),
                        threshold=0.7,
                        orgdb="org.Hs.eg.db")
      
      p <-
        plot_mat2 %>%
        as_tibble() %>%
        group_by(parentTerm) %>%
        mutate(n = row_number()) %>%
        ggplot(aes(area = size, subgroup = parentTerm)) +
        
        geom_treemap(aes(fill = parentTerm,
                         alpha = n),
                     color = "black",
                     show.legend = F) +
        geom_treemap_subgroup_border(color = "black") +
        
        geom_treemap_text(aes(label = term),
                          colour = "black",
                          place = "centre",
                          alpha = 0.4,
                          grow = TRUE) +
        
        geom_treemap_subgroup_text(place = "centre",
                                   grow = T,
                                   reflow = T,
                                   alpha = 1,
                                   colour = "white",
                                   fontface = "bold",
                                   min.size = 0) +
        scale_alpha_continuous(range = c(0.8, 1)) +
        scale_fill_manual(values = rep(ggthemes::tableau_color_pal(palette = "Hue Circle",
                                                                   type = c("regular"),
                                                                   direction = 1)(19),
                                       10)) +
        theme(plot.margin = c(1.5, 1.5, 1.5, 1.5)) +
        theme_void() +
        ggtitle(GO_ontology) + 
        theme(plot.margin = margin(5.5, 5.5, 5.5, 10))
    } else {
      treemap <- 
        ggplot() +
        theme_void() +
        ggtitle(GO_ontology)
    }
    
  }




annotation_report <- function(i) {
  
  genes <- final_clustering %>% 
    filter(value==i) %>% 
    pull(gene)
  
  n_genes <- length(genes)
  
  p_h <- 
    scaled_data %>% 
    filter(enssscg_id %in% genes) %>% 
    select(enssscg_id,tissue,scaled_value) %>%
    spread(key=tissue, value = scaled_value) %>%
    column_to_rownames("enssscg_id") %>% 
    pheatmap(clustering_method = "ward.D2", 
             color = viridis(20, direction = -1, option = "B"),
             show_rownames = F,border_color = NA, fontsize = 6)
  
  p_bp <-
    fisher_res %>% 
    filter(cluster ==i) %>% 
    filter(adj_p < 0.05) %>% 
    ggplot(aes(x=odds_ratio, y = enhanced_tissues,  fill=enhanced_tissues)) +
    geom_bar(stat = "identity", show.legend = F, width = 0.5) +
    theme_minimal() +
    scale_fill_manual(values = tissue_colors)
 
   p_d <- data_tissue %>% 
    filter(enssscg_id %in% genes) %>% 
    ggplot(aes(x=log10(tissue_exp + 1))) +
    geom_density() +
    theme_minimal()
  
  p_D <- 
    scaled_data %>% 
    filter(enssscg_id %in% genes) %>% 
    ggplot(aes(x=log10(scaled_value)))
  
  p_tm_bp <- plot_treemap(GO_terms = GO_enrichments[[1]], GO_ontology = "BP", clus = i)
  p_tm_cc <- plot_treemap(GO_terms = GO_enrichments[[2]], GO_ontology = "CC", clus = i)
  p_tm_mf <- plot_treemap(GO_terms = GO_enrichments[[3]], GO_ontology = "MF", clus = i)
  
  # 
  # as.ggplot(p_h) / ((p_bp / p_d) | (p_tm_bp / p_tm_cc / p_tm_mf)) +
  #                     plot_annotation(title = paste("Cluster",i, ":  n =", n_genes))
  # 
  as.ggplot(p_h) / ((p_bp | p_d) / (p_tm_bp | p_tm_cc | p_tm_mf)) +
    plot_annotation(title = paste("Cluster",i, ":  n =", n_genes)) 
}

