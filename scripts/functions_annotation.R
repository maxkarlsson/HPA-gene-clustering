# Functopm for hypergeometric test

gene_fischer_test <-
  function(cluster_membership, term_membership, alternative = "greater") {
    
    # cluster_membership <- clust_list 
    # term_membership <- tissue_list 
    n_genes <-
      term_membership %>%
      unlist() %>%
      n_distinct()
    
    fisher_data <-
      cluster_membership %>%
      enframe("gene", "cluster") %>%
      full_join(term_membership %>%
                  map(enframe, "row", "gene") %>%
                  bind_rows(.id = "term")) %>%
      mutate(term = factor(term),
             cluster = factor(cluster)) %>%
      group_by(cluster, term, .drop = F) %>%
      count() %>%
      group_by(term) %>%
      mutate(total_term = sum(n),
             total_nonterm = n_genes - total_term) %>%
      filter(!is.na(cluster) &
               !is.na(term)) %>%
      group_by(cluster) %>%
      mutate(total_in_cluster = sum(n),
             total_not_in_cluster = n_genes - total_in_cluster ) %>%
      ungroup()
    
    
    
    fisher_data %>%
      group_by(cluster, term) %>%
      do({
        g_data <<- .
        
        contingency_matrix <-
          matrix(c(g_data$n,
                   g_data$total_in_cluster - g_data$n,
                   g_data$total_term - g_data$n,  
                   g_data$total_nonterm - (g_data$total_term - g_data$n)),
                 nrow = 2,
                 dimnames =
                   list(c("Term", "Nonterm"),
                        c("In cluster", "Not in cluster")))
        
        contingency_matrix %>%
          fisher.test(alternative = alternative) %>%
          tidy() %>%
          mutate(term_overlap = paste(contingency_matrix[1],
                                      contingency_matrix[1] + contingency_matrix[3],
                                      sep = "/"),
                 cluster_membership = paste(contingency_matrix[1],
                                            contingency_matrix[1] + contingency_matrix[2],
                                            sep = "/"))
      }) %>%
      ungroup() %>%
      mutate(p.value = ifelse(p.value == 0, .Machine$double.xmin, p.value),
             adj_p = p.adjust(p.value, method = "BH")) %>%
      rename(odds_ratio = estimate) %>%
      arrange(adj_p)
  }



cluster_query <-
  function(query_data, clustering, distance_method = "euclidean") {
    
    stopifnot(clustering %in% c("ward.D", "ward.D2", "single", "complete",
                                "average", "mcquitty", "median", "centroid"))
    
    query_data %>%
      select(gene_name, tissue, exp) %>%
      mutate(exp = ifelse(!is.finite(exp), 0, exp)) %>%
      spread(tissue, exp) %>%
      column_to_rownames("gene_name") %>%
      
      {list(gene = dist(., distance_method) %>%
              stats::hclust(method = clustering),
            sample = t(.) %>%
              dist(distance_method) %>%
              stats::hclust(method = clustering))}
    
  }

# Function to generate treemap

plot_treemap <- 
  function(GO_terms,GO_ontology, clus) {
    plot_data <-
      GO_terms$enriched_terms %>% 
      filter(Cluster == clus) %>%
      filter(qvalue < 0.05)#0.05
    
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


# Function to generate annotaiton report

annotation_report <- function(i) {
  
  genes <- final_clustering %>% 
    filter(value==i) %>% 
    pull(gene)
  
  n_genes <- length(genes)
  
  
  gene_ann <- 
    row_annotation %>% 
    filter(ensg_id %in% genes)
  
  p_h <- 
    scaled_data %>% 
    filter(enssscg_id %in% genes) %>% 
    select(enssscg_id,tissue,scaled_value) %>%
    spread(key=tissue, value = scaled_value) %>%
    column_to_rownames("enssscg_id") %>% 
    pheatmap(clustering_method = "ward.D2", 
             color = viridis(20, direction = -1, option = "B"),
             show_rownames = F,border_color = NA, fontsize = 5,
             annotation_row = gene_ann  %>%  column_to_rownames("ensg_id"), 
             annotation_col = tissue_ann,
             annotation_colors = all_ann, annotation_legend = F)
  
  p_h_bp <- 
    scaled_data %>% 
    filter(enssscg_id %in% genes) %>% 
    group_by(tissue) %>% 
    mutate(sumed_expression = mean(scaled_value)) %>% 
    ungroup() %>%
    select(tissue,sumed_expression) %>% 
    left_join(colors %>%  select(tissue,tissue_name, consensus_tissue_name, organ_name)) %>% 
    distinct() %>% 
    arrange(organ_name) %>%
    mutate(tissue_name = factor(tissue_name, unique(tissue_name))) %>%
    ggplot(aes(x= tissue_name, y = sumed_expression, fill = consensus_tissue_name)) +
    geom_bar(stat = "identity", show.legend = F) +
    scale_fill_manual(values = tissue_colors) +
    scale_y_continuous(expand = expansion(c(0, 0.05)), limits = c(0,1)) +
    theme_minimal() +
    # scale_x_discrete(breaks = NULL)
    theme(panel.grid.major.x = element_blank() ,
          panel.grid.major.y = element_line( size=.1 ),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
          axis.text.y = element_text(size = 6)) +
    ylab("Average relative expression") +
    xlab("") +
    ggtitle("Average tissue expression") + 
    theme(plot.title = element_text(size=10))
  
  enriched_tissues <- 
    fisher_res %>% 
    filter(cluster == i) %>% 
    filter(adj_p < 0.05) %>% 
    pull(enhanced_tissues) %>% 
    as.character()
  
  if(length(enriched_tissues) > 0) {
    percentage_enrichment <- 
      tissue_info %>% 
      filter(ensg_id %in% genes) %>% 
      filter(enhanced_tissues %in% enriched_tissues) %>% 
      group_by(enhanced_tissues) %>% 
      count() %>% 
      mutate(n = n/length(genes))
    
    p_bp <-
      fisher_res %>% 
      filter(cluster ==i) %>% 
      filter(adj_p < 0.05) %>% 
      left_join(percentage_enrichment) %>% 
      ggplot(aes(x=odds_ratio, y = enhanced_tissues,  fill=enhanced_tissues)) +
      geom_bar(stat = "identity", show.legend = F, width = 0.5) +
      theme_minimal() +
      scale_fill_manual(values = tissue_colors) +
      theme(plot.title = element_text(size=10),
            axis.text = element_text(size = 6)) +
      ggtitle("Enrichment towards tissuespecificity annotation") +
      geom_text(aes(label = paste(round(n,2)*100, "%")), position = position_dodge(width = 0.8), vjust = -0.6) 
  }
  
  else {
    p_bp <- 
      ggplot() +
      theme_void() +
      ggtitle("Enrichment towards tissuespecificity annotation") +
      theme(plot.title = element_text(size=10),
            axis.text = element_text(size = 6))
    
  }
  
  
  p_d <- data_tissue %>% 
    filter(enssscg_id %in% genes) %>% 
    ggplot(aes(x=log10(tissue_exp + 1))) +
    geom_density() +
    theme_minimal() +
    ggtitle("Density plot") + 
    theme(plot.title = element_text(size=10),
          axis.text = element_text(size = 6))
  
  
  p_tm_bp <- plot_treemap(GO_terms = GO_enrichments[[1]], GO_ontology = "BP", clus = i)
  p_tm_cc <- plot_treemap(GO_terms = GO_enrichments[[2]], GO_ontology = "CC", clus = i)
  p_tm_mf <- plot_treemap(GO_terms = GO_enrichments[[3]], GO_ontology = "MF", clus = i)
  
  # 
  # as.ggplot(p_h) / ((p_bp / p_d) | (p_tm_bp / p_tm_cc / p_tm_mf)) +
  #                     plot_annotation(title = paste("Cluster",i, ":  n =", n_genes))
  # 
  as.ggplot(p_h) / (p_h_bp /(p_bp | p_d) / (p_tm_bp | p_tm_cc | p_tm_mf)) +
    plot_annotation(title = paste("Cluster",i, ":  n =", n_genes)) +
    # plot_layout(guides = 'collect')
    plot_layout(widths = c(3, 2,2), heights = c(2,4))
}

