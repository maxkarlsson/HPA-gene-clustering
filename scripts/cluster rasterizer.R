

library(tidyverse)
library(magrittr)
library(patchwork)
library(pheatmap)
library(sf)
library(sp)
library(concaveman)
library(rayshader)


get_density <- 
  function(x, y, h = 0.5, n = 100, lims = c(range(x),
                                            range(y))) {
    g_density <- 
      MASS::kde2d(x, y, h = h, n = n, 
                  lims = lims) 
    
    
    g_density$z %>%
      as_tibble(rownames = "x") %>%
      gather(y, z, -x) %>% 
      mutate(y = gsub("V", "", y)) %>% 
      mutate_at(vars(x,y), as.integer) %>%
      left_join(enframe(g_density$x, "x", "x_coord"),
                by = "x") %>% 
      left_join(enframe(g_density$y, "y", "y_coord"),
                by = "y")
  }


generate_cluster_hulls <- 
  function(V1,
           V2, 
           element_id,
           cluster_membership, 
           n = 300,
           cum_z_lim = 0.95,
           frac_lim = 0.05,
           plot_range = c(range(V1),
                          range(V2)) * 1.05,
           poly_concavity = 1,
           poly_smoothing = 1,
           relative_bandwidth = 1/200) {
    require(tidyverse)
    require(magrittr)
    require(sf)
    require(sp)
    require(concaveman)
    
    # Compile indata
    cluster_data <- 
      tibble(V1, V2, 
             element_id,
             cluster = cluster_membership)
    
    # Calculate plot range
    plot_range_tb <- 
      set_names(plot_range,
                c("xmin", 
                  "xmax",
                  "ymin", 
                  "ymax")) %>% 
      enframe() %>% 
      spread(name, value)
    
    # Calculate diagonal length
    plot_diagonal <- 
      sqrt((plot_range_tb$xmax - plot_range_tb$xmin)^2 + 
             (plot_range_tb$ymax - plot_range_tb$ymin)^2)
    
    # Set bandwidth to fraction of diagonal
    plot_bandwidth <- 
      relative_bandwidth * plot_diagonal
    
    # Find subclusters
    subclusters <-
      cluster_data %>%
      group_by(cluster) %>%
      mutate(n_cluster_genes = n_distinct(element_id), 
             sub_cluster = data.frame(V1, V2) %>%
               fpc::dbscan(eps = plot_bandwidth) %$%
               cluster) %>%
      ungroup() %>% 
      group_by(cluster, sub_cluster) %>%
      mutate(n_sub_genes = n_distinct(element_id)) %>% 
      ungroup() 
    
    # Classify subclusters
    subclusters_classes <- 
      subclusters %>% 
      select(cluster, sub_cluster, n_cluster_genes, n_sub_genes) %>% 
      distinct() %>% 
      group_by(cluster) %>%
      mutate(sub_type = case_when(sub_cluster == 0 ~ "outlier",
                                  n_sub_genes / n_cluster_genes < frac_lim ~ "outlier",
                                  rank(-n_sub_genes, 
                                       ties.method = "first") == 1 ~ "primary",
                                  T ~ "secondary")) %>% 
      select(cluster, sub_cluster, sub_type)
    
    subclusters_classed <- 
      subclusters %>% 
      left_join(subclusters_classes,
                by = c("cluster", "sub_cluster"))
    
    
    # Calculate plot density
    plot_density <- 
      subclusters_classed %>% 
      filter(sub_type != "outlier") %>%
      group_by(cluster, sub_cluster) %>%
      do({
        get_density(.$V1, 
                    .$V2, 
                    h = plot_bandwidth, 
                    n = n, 
                    lims = plot_range) 
        
      }) %>% 
      ungroup() %>% 
      filter(z > 1e-200) %>% 
      group_by(cluster, sub_cluster) %>% 
      mutate(z = z / sum(z)) %>% 
      arrange(cluster, sub_cluster, -z) %>% 
      mutate(cum_z = cumsum(z)) %>% 
      ungroup()
    
    
    # Filter pixels such that 95% of density is included
    # Each point is then assigned to the cluster with highest density
    plot_density_filtered <- 
      plot_density %>% 
      filter(cum_z < cum_z_lim) %>%
      group_by(x, y) %>% 
      top_n(1, z) %>%
      slice(1) %>% 
      ungroup()
    
    
    # Calculate size of landmass
    plot_density_landmass <-
      plot_density_filtered %>%
      group_by(cluster, sub_cluster) %>%
      mutate(landmass = data.frame(x_coord, y_coord) %>%
               fpc::dbscan(eps = plot_bandwidth) %$%
               cluster) %>%
      group_by(cluster, sub_cluster, landmass) %>% 
      mutate(n_landmass_points = length(x)) %>%  
      ungroup() %>% 
      group_by(cluster) %>%
      mutate(n_total_points = length(x)) %>% 
      ungroup() 
    
    # Classify landmasses
    plot_density_landmass_classes <- 
      plot_density_landmass %>% 
      select(cluster, sub_cluster, landmass, n_landmass_points, n_total_points) %>% 
      distinct() %>%
      group_by(cluster, sub_cluster) %>%
      
      mutate(frac_landmass = n_landmass_points / n_total_points,
             landmass_type = case_when(rank(-n_landmass_points, 
                                            ties.method = "first") == 1 ~ "primary",
                                       T ~ "secondary")) %>% 
      ungroup() %>% 
      select(cluster, sub_cluster, landmass, landmass_type, frac_landmass)
    
    plot_density_landmass_classed <- 
      plot_density_landmass %>%
      left_join(plot_density_landmass_classes,
                by = c("cluster", 
                       "sub_cluster",
                       "landmass")) %>% 
      arrange(cluster)
    
    plot_density_mainland_filtered <-
      plot_density_landmass_classed %>% 
      filter(frac_landmass > frac_lim) %>%
      filter(cum_z < cum_z_lim) %>%
      group_by(x, y) %>% 
      top_n(1, z) %>%
      slice(1) %>% 
      ungroup()
    
    ######
    
    # Create polygons
    # Poly smoothing: How small distances should be further detailed - 
    # higher values --> less detailed
    # poly concavity: How convex polygons should be - 
    # higher values --> less detailed
    
    plot_data_hulls <- 
      plot_density_mainland_filtered %>% 
      
      group_by(cluster, sub_cluster, landmass) %>% 
      do({
        st_as_sf(., coords=c('x_coord','y_coord')) %>%
          concaveman(concavity = poly_concavity, 
                     length_threshold = plot_bandwidth * poly_smoothing) %$%
          st_coordinates(polygons) %>% 
          as_tibble()
      }) %>% 
      ungroup() %>% 
      mutate(polygon_id = paste(cluster, sub_cluster, landmass, sep = "_"))
    
    
    
    list(hulls = plot_data_hulls,
         density = plot_density,
         landmass_pixels = plot_density_mainland_filtered)
  }



