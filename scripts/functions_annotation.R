# Functopm for hypergeometric test

gene_fischer_test <-
  function(cluster_membership, term_membership, alternative = "greater") {
    require(broom)
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



# Formatting functions:


GO_formatting_function <- 
  . %>% 
  select(ensg_id = `Gene stable ID`,
         term = `GO term name`,
         term_id = `GO term accession`) %>% 
  filter(!is.na(term)) %>% 
  distinct() %>% 
  filter(ensg_id %in% clustering_data$gene) %>% 
  group_by(term_id, term) %>% 
  filter(n_distinct(ensg_id) >= 8) %>% 
  ungroup()

HPAspecificity_formatting_function <- 
  . %>% 
  select(ensg_id,
         term = enhanced_tissues) %>% 
  mutate(term = gsub(", ", ";", term)) %>% 
  separate_rows(term, sep = ",") %>% 
  mutate(term = gsub(";", ", ", term)) %>% 
  filter(!is.na(term)) %>% 
  filter(term != "NULL") %>% 
  mutate(term = trimws(term), 
         term_id = term)


database_formatting_functions <- 
  list(secretome_location = . %>% 
         select(ensg_id = Ensembl,
                term = `Secretome location`) %>% 
         filter(!is.na(term)) %>% 
         mutate(term_id = term),
       
       specificity_blood = HPAspecificity_formatting_function,
       specificity_brain = HPAspecificity_formatting_function,
       specificity_tissue = HPAspecificity_formatting_function,
       specificity_celline = HPAspecificity_formatting_function,
       specificity_singlecell = HPAspecificity_formatting_function,
       
       subcellular_location = . %>% 
         filter(Gene %in% clustering_data$gene) %>% 
         select(ensg_id = Gene, 
                Enhanced, 
                Supported, 
                Approved) %>% 
         gather(reliability, location, -1) %>% 
         filter(!is.na(location)) %>% 
         separate_rows(location, sep = ";") %>% 
         select(ensg_id, term = location) %>% 
         mutate(term_id = term),
       
       panglao_cellmarkers = . %>% 
         filter(grepl("Hs", species), 
                `gene type` %in% c("protein-coding gene",
                                   "protein coding gene")) %>% 
         select(gene_name = 2, 
                term = 3) %>% 
         distinct() %>% 
         mutate(term_id = term),
       
       protein_class = . %>% select(ensg_id = Ensembl,
                                    protein_class = `Protein class`) %>% 
         separate_rows(protein_class, sep = ", ") %>% 
         rename(term = protein_class) %>% 
         mutate(term_id = term),
       
       reactome = . %>% 
         filter(species == "Homo sapiens") %>% 
         select(ensg_id, term = description, term_id = id) %>% 
         distinct() %>% 
         filter(ensg_id %in% clustering_data$gene), 
       
       trrust = . %>% 
         select(gene_name, term) %>% 
         mutate(term_id = term),
       
       GO_BP_original = . %>% 
         filter(`GO domain` == "biological_process") %>% 
         GO_formatting_function,
       GO_MF_original = . %>% 
         filter(`GO domain` == "molecular_function") %>% 
         GO_formatting_function,
       GO_CC_original = . %>% 
         filter(`GO domain` == "cellular_component") %>% 
         GO_formatting_function)



gene_formatting_functions <- 
  list(`gene names` = . %>% 
         left_join(gene_info %>% 
                     select(ensg_id, gene_name),
                   by = c("gene_name")) %>% 
         select(ensg_id, term, term_id))


get_db <- 
  function(enrichment_settings, db_id) {
    
    db_file <- 
      enrichment_settings %>% 
      filter(id == db_id) %>% 
      pull(file)
    
    raw_db <- 
      read_tsv(paste0("annotation_databases/", db_file))
    
    db_format_function <- 
      database_formatting_functions[[db_id]]
    
    raw_db %>% 
      db_format_function()
  }

get_specificity_db <- 
  function(enrichment_settings, db_id) {
    
    db_file <- 
      enrichment_settings %>% 
      filter(id == db_id) %>% 
      pull(file)
    
    raw_db <- 
      read_tsv(paste0("annotation_databases/", db_file))
    
    
    out_db <- 
      raw_db %>% 
      set_colnames(c("ensg_id",
                     "specificity_category",
                     "distribution_category",
                     "ts_score",
                     "enhanced_score",
                     "enhanced_tissues"))
    
    out_db
  }


perform_ORA <-
  function(gene_associations,
           database,
           universe,
           n_clusters = 5,
           minGSSize = 10,
           maxGSSize = Inf,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2) {
    require(clusterProfiler)
    require(multidplyr)
    
    if(n_clusters != 1) {
      worker_cluster <- new_cluster(n = n_clusters)
      cluster_library(worker_cluster, c("dplyr",
                                        "tidyverse"))
      cluster_copy(worker_cluster, c("enricher",
                                     "universe",
                                     "database",
                                     "minGSSize",
                                     "maxGSSize",
                                     "pvalueCutoff",
                                     "qvalueCutoff" ))
      
      pre_out <- 
        gene_associations %>%
        group_by(partition) %>%
        partition(worker_cluster) 
    } else {
      pre_out <- 
        gene_associations %>%
        group_by(partition)
    }
    
    outdata <-
      pre_out %>% 
      do({
        g_data <- .
        pull(g_data, gene) %>%
          enricher(universe = universe,
                   TERM2GENE = database, 
                   minGSSize = minGSSize,
                   maxGSSize = maxGSSize,
                   pvalueCutoff = pvalueCutoff,
                   qvalueCutoff = qvalueCutoff) %>%
          as_tibble()
      }) %>%
      ungroup() %>%
      collect()
    
    if(n_clusters != 1) rm(worker_cluster)
    outdata
  }
