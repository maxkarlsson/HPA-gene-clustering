
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


# Plot PCA results

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


# UMAP function

umap_calc <- function (data, row_names, n_neigh, n_comp = 2, n_ep = 1000, color_tissue = F, met = "euclidean") {
  set.seed(42)
  rows <- data[[row_names]]
  umap_res <-
    data %>%
    column_to_rownames(row_names) %>%
    uwot::umap(n_neighbors = n_neigh, n_components = n_comp, n_epochs = n_ep, nn_method = "annoy", metric = met) %>%
    tibble::as_tibble () %>%
    mutate(features = rows)
  
  return(umap_res)
}


get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


# Plot UMAP results

umap_plot <- function (umap_res, color_by = NULL, color_groups = NULL) {
  if (missing(color_by)) {
    uplot <-
      umap_res %>%
      ggplot(aes(V1,V2)) +
      geom_point() +
      theme_minimal()+
      coord_fixed()
  }
  
  else if (color_by == "density") {
    umap_res$density <- get_density(umap_res$V1, umap_res$V2, n = 1000)
    uplot <-
      umap_res %>%
      ggplot(aes(V1,V2, colour = density)) +
      geom_point() +
      theme_minimal()+
      coord_fixed()+ 
      scale_color_viridis() 
    }
  
  else if (color_by == "tissue") {
    uplot <-
      umap_res %>%
      mutate(enssscg_id = gene_names) %>%
      left_join(color_groups, by = c("enssscg_id" = "ensg_id")) %>%
      ggplot(aes(V1,V2, color = enhanced_tissues)) +
      geom_point(show.legend = F) +
      theme_minimal() +
      scale_color_manual(values = tissue_colors) +
      coord_fixed()
    }
  
  else if (color_by == "cluster") {
    num_clusters <- length(unique(color_groups$value))
    pal <- colorRampPalette(tableau_color_pal(palette = "Classic Cyclic")(13))(num_clusters)
    uplot <-
      umap_res %>%
    #  mutate(enssscg_id = gene_names) %>%
      left_join(color_groups, by = c("features" = "gene")) %>%
      ggplot(aes(V1,V2, color = as.factor(value))) +
      geom_point(show.legend = F, size = 0.3) +
      theme_simple +
      scale_color_manual(values = pal) +
      coord_fixed()
    
    }
  
    else {
    uplot <-
      umap_res %>%
      ggplot(aes(V1,V2)) +
      geom_point() +
      theme_minimal()+
      coord_fixed()
  }

    return(uplot)
}


# Plot clusteirng solution on UMAP

clustering_umaps <- function(name, umaps, partitions) { #, 
  partition <-
    partitions %>% 
    filter(id == name)
  
  if(grepl("zscore euclidean",name)) {
    umap <- umaps$`zscore euclidean`
  }
  
  else if(grepl("zscore pearson",name)) {
    umap <- umaps$`zscore correlation`
  }
  else if(grepl("max euclidean",name)) {
    umap <- umaps$`max euclidean`
  }
  else if(grepl("max pearson",name)) {
    umap <- umaps$`max correlation`
  }
  else if(grepl("random",name)) {
    umap <- umaps$`zscore correlation`
  }
  
  
  umap_plot(umap_res = umap, color_by = "cluster", color_groups = partition) +
    ggtitle(name)
}
