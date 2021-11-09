# ----- 3D UMAP -----
graph_umap_3d <- 
  data_cluster$neighbors$nn %>%
  {. + t(.)} %>%
  umap(n_components = 3)

umap_3d %>% 
  as_tibble() %>% 
  mutate(gene = rownames(data_cluster$pca_data$scores)) %>% 
  left_join(data_cluster$cluster_data %>% 
              filter(resolution == 2, seed == 1)) %>%
  plot_ly(x = ~V1, 
          y = ~V2,
          z = ~V3,
          color = ~cluster,
          colors = cluster_palette, 
          size = 0.1)

# umap_3d %>% 
#   as_tibble() %>% 
#   mutate(gene = rownames(data_cluster$pca_data$scores)) %>% 
#   left_join(data_cluster$cluster_data %>% 
#               filter(resolution == 2, seed == 1)) %>%
#   plot_ly(x = ~V1, 
#           y = ~V2,
#           z = ~V3,
#           color = ~gene,
#           colors = set_names(ensg_exp_palette$color,
#                              ensg_exp_palette$ensg_id), 
#           size = 0.1,
#           showscale=FALSE) %>% 
#   add_trace(showlegend = FALSE)
#   add_trace(showscale = FALSE)

# graph_umap_3d <- 
#     data_cluster$neighbors$nn %>%
#     {. + t(.)} %>%
#     umap(n_components = 3, n_neighbors = 15)

p <- 
  graph_umap_3d %>% 
  as_tibble() %>% 
  mutate(gene = rownames(data_cluster$pca_data$scores)) %>% 
  left_join(data_cluster$cluster_data %>% 
              filter(resolution == 2, seed == 1)) %>%
  mutate(cluster = factor(cluster)) %>%
  plot_ly(x = ~V1, 
          y = ~V2,
          z = ~V3,
          color = ~cluster,
          colors = set_names(cluster_palette,
                             0:(length(cluster_palette) - 1)), 
          size = 0.1)

p

saveWidgetFix(p, "results/3d_umap.html", selfcontained = T, libdir = "lib")

tissue_contr_data <-
  data %>% 
  group_by(ensg_id) %>% 
  mutate(tmm = tmm / sum(tmm))



gene_info92 <- read_tsv("data/meta/geninfo_92.tsv")

data %>%
  inner_join(data_cluster$cluster_data %>% 
               filter(resolution == 2, seed == 1) %>% 
               filter(cluster == 43),
             by = c("ensg_id" = "gene")) %>%
  left_join(gene_info92) %>%
  select(celltype, tmm, gene_name) %>%
  spread(celltype, tmm) %>% 
  column_to_rownames("gene_name") %>% 
  {log10(. + 1)} %>%
  pheatmap()

data_cluster$cluster_data %>% 
  filter(resolution == 2, seed == 1) %>% 
  filter(cluster == 43) %>%
  left_join(gene_info92,
            by = c("gene" = "ensg_id")) %>%
  pull(gene_name) %>% 
  paste(collapse = "\n") %>% 
  cat()

tissue_contr_data$celltype %>% unique
graph_umap_3d %>% 
  as_tibble() %>% 
  mutate(gene = rownames(data_cluster$pca_data$scores)) %>% 
  left_join(data_cluster$cluster_data %>% 
              filter(resolution == 0.6, seed == 1)) %>%
  left_join(tissue_contr_data %>% 
              filter(celltype == "testis"),
            by = c("gene" = "ensg_id")) %>%
  plot_ly(x = ~V1, 
          y = ~V2,
          z = ~V3,
          name = ~cluster,
          color = ~tmm,
          # colors = cluster_palette, 
          size = 0.1)