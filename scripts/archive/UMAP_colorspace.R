
shared_colorspace_facet_plot <- 
  function(data, V1, V2, point_id, facet_id, ref_id) {
    plot_data <- 
      data %>% 
      rename(V1 = V1, 
             V2 = V2,
             point_id = point_id, 
             facet_id = facet_id) 
    
    
    ref_palette <- 
      plot_data %>% 
      filter(facet_id == ref_id) %>% 
      mutate(r = scales::rescale(V1, c(0,1)),
             g = scales::rescale(V2, c(0,1)),
             b = 0.5,
             color = rgb(r, g, b)) %>% 
      select(point_id, color)
    
    plot_data %>% 
      left_join(ref_palette,
                by = "point_id") %>% 
      mutate(color = ifelse(is.na(color), 
                            "darkgray", 
                            color))
    
  }

plot_data <- 
  all_graph_umaps %>% 
  select(dataset_id) %>% 
  distinct() %>% 
  group_by(dataset_id) %>% 
  do({
    shared_colorspace_facet_plot(data = all_graph_umap_data_scaled,
                                 V1 = "UMAP_1",
                                 V2 = "UMAP_2",
                                 point_id = "gene",
                                 facet_id = "dataset_id",
                                 ref_id = .$dataset_id)
  }) %>% 
  ungroup()


plot_data %>% 
  ggplot(aes(V1, V2, color = color)) +
  geom_point(size = 0.1) +
  theme_bw() +
  coord_fixed() +
  facet_grid(dataset_id~facet_id) +
  scale_color_identity()


plot_data %>% 
  group_by(dataset_id, 
           facet_id) %>% 
  
  do({
    g_data <- .
    g_data %$% 
      makeHexData(x = V1, 
                  y = V2, 
                  z = color, 
                  bins = 50, color_mean)
  }) %>% 
  ggplot(aes(x = x, y = y, fill = z, group = paste(dataset_id, facet_id))) +
  geom_hex(stat = "identity") +
  facet_grid(dataset_id~facet_id) +
  theme_bw() +
  coord_fixed() +
  scale_fill_identity()
