
# PCA function

pca_calc <- function(data, npcs, m = NULL) {
  require(pcaMethods)
  pca_res <-
    data %>%
    pca(nPcs = npcs)
  
  pca_stats <-
    summary(pca_res) %>%
    t() %>%
    as_tibble(rownames = "PC") %>%
    mutate(PC = as.numeric(gsub("PC", "", PC))) %>%
    dplyr::rename(R2cum = 3)
  
  informative_pcs <- pca_stats$PC[which(pca_stats$R2cum > 0.95)[1]]
  
  pca_stats_plot <-
    pca_stats %>%
    select(PC, R2cum) %>%
    ggplot(aes(PC, R2cum)) +
    geom_point() +
    geom_line() +
    theme_minimal() +
    geom_vline(xintercept = informative_pcs, linetype = "dashed") +
    annotate("text",
             x = informative_pcs,
             y = 0.55,
             label = paste0("PC ", informative_pcs,
                            "\nR2 = ", round(pca_stats[informative_pcs, ]$R2cum, 3)),
             hjust = 1,
             vjust = 0)
  
  
  if(!is.null(m)) {
    return(list(pca = pca_res,
       scores = pcaMethods::scores(pca_res),
       loadings = pcaMethods::loadings(pca_res),
       stats = pca_stats,
       pc_95 = informative_pcs,
       plot = pca_stats_plot,
       method = m))
    }
  
  else {
    return(list(pca = pca_res,
                scores = pcaMethods::scores(pca_res),
                loadings = pcaMethods::loadings(pca_res),
                stats = pca_stats,
                pc_95 = informative_pcs,
                plot = pca_stats_plot))
    }
  }

# Plot PCA

pca_plot <- function(data) {
  data$scores %>%
    as_tibble()%>%
    mutate(enssscg_id = scaled_data$enssscg_id) %>%
    left_join(tissue_info, by = c("enssscg_id" = "ensg_id")) %>%
    ggplot(aes(PC1,PC2, color = enhanced_tissues)) +
    geom_point(show.legend = F) +
    theme_minimal()+
    scale_color_manual(values = tissue_colors)
  
}





# UMAP function (?)

umap_calc <- function (data, row_names, n_neigh, n_comp = 2, n_ep = 1000, color_tissue = F) {
  set.seed(42)
  rows <- data[[row_names]]
  umap_res <-
    data %>%
    column_to_rownames(row_names) %>%
    uwot::umap(n_neighbors = n_neigh, n_components = n_comp, n_epochs = n_ep) %>%
    tibble::as_tibble () %>%
    mutate(features = rows)
  
  
  if (color_tissue == T) {
    umap_plot <-
      umap_res %>%
      mutate(enssscg_id = gene_names) %>%
      left_join(tissue_info, by = c("enssscg_id" = "ensg_id")) %>%
      ggplot(aes(V1,V2, color = enhanced_tissues)) +
      geom_point(show.legend = F) +
      theme_minimal() +
      scale_color_manual(values = tissue_colors) +
      coord_fixed()
  }
  
  else {
    umap_plot <-
      umap_res %>%
      ggplot(aes(V1,V2)) +
      geom_point() +
      theme_minimal()+
      coord_fixed()
  }
  
  list(umap = umap_res,
       plot = umap_plot)
}
