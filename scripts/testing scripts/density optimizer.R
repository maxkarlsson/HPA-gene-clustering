
plot_cluster <- 11

plot_density %>% 
  filter(cluster == plot_cluster) %>% 
  ggplot(aes(x_coord, y_coord, fill = z)) +
  geom_tile() +
  geom_point(data = subclusters_classed %>% 
               filter(cluster == plot_cluster),
             aes(V1, V2, color = sub_type), 
             inherit.aes = F)
  


pixel_dist_data <- 
  subclusters_classed %>% 
  filter(sub_type != "outlier") %>% 
  filter(cluster == plot_cluster) %>% 
  select(V1, V2, element_id, cluster, sub_cluster) %>% 
  left_join(plot_density %>% 
              select(cluster, sub_cluster, x_coord, y_coord, z)) %>% 
  mutate(dist = sqrt((V1 - x_coord) ^ 2 + (V2 - y_coord) ^ 2)) %>% 
  group_by(x_coord, y_coord, z) %>% 
  summarise(min_dist = min(dist)) %>% 
  ungroup()
  
pixel_dist_data %>% 
  arrange(min_dist) %>% 
  mutate(cumz = cumsum(z), 
         bws = min_dist/plot_bandwidth) %>% 
  # View
  ggplot(aes(bws, cumz)) +
  geom_hex(bins = 100)


plot_lims <- 
  subclusters_classed %>% 
  filter(sub_type != "outlier") %>% 
  filter(cluster == plot_cluster) %>% 
  summarise(xmin = min(V1) - plot_bandwidth,
            xmax = max(V1) + plot_bandwidth,
            ymin = min(V2) - plot_bandwidth,
            ymax = max(V2) + plot_bandwidth)

plot_density %>% 
  filter(cluster == plot_cluster) %>% 
  filter(x_coord > plot_lims$xmin,
         x_coord < plot_lims$xmax,
         y_coord > plot_lims$ymin,
         y_coord < plot_lims$ymax) %>% 
  pull(z) %>% 
  sum


un_pixs <-
  plot_density %>% 
  select(x_coord, y_coord) %>% 
  distinct() 

un_pixs %>% 
  # filter(!y_coord %in% seq(plot_range[3], plot_range[4], length.out = n))
  filter(!x_coord %in% seq(plot_range[1], plot_range[2], length.out = n))

pixels_x <- 
  seq(plot_range[1], plot_range[2], length.out = n)

pixels_y <- 
  seq(plot_range[3], plot_range[4], length.out = n)


x_step <- pixels_x[1] - pixels_x[2]
y_step <- pixels_y[1] - pixels_y[2]



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
