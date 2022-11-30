
# ------ Load stuff ----

library(tidyverse)
library(magrittr)
library(patchwork)
library(ggrastr)
library(ggplotify)
library(readxl)

source("scripts/functions_utility.R")

region_colors <- 
  read_tsv("data/colors/colors_regions.tsv")

region_palette <- 
  region_colors %>% 
  select(-1) %>% 
  distinct() %>% 
  deframe() %>% 
  c(c("female organs" = .[["vagina"]],
      "pancreatic cells" = .[["pancreas"]],
      "monocytes & neutrophils" = .[["NK-cells"]],
      "nk-cells & t-cells" = .[["NK-cells"]],
      "ciliated cells" = .[["Basal respiratory cells"]],
      "spermatogonia & spermaocytes" = .[["testis"]],
      "spermaocytes & spermatogonia" = .[["testis"]],
      "spermaocytes & erythroid cells" = .[["testis"]],
      "Plasma cells (B-cell immune response cascade)" = .[["NK-cells"]],
      "epithelial cell types" = .[["Squamous epithelial cells"]],
      "plasmacytoid DCs" = .[["NK-cells"]]))

region_palette[["plasmacytoid DCs"]] <- 
  region_palette[["NK-cells"]]

region_palette
# region_palette[["female organs"]] %>% prismatic::color()

consensus_colors <- 
  read_tsv("data/colors/colors_consensus.tsv")

###

dataset_metadata <- read_csv("run_settings/20221004_settings.csv") %>% 
  filter(dataset != "brain")


data_paths <- read_csv("run_settings/20221020 all_datasets.csv")

enrichment_settings <- read_csv("run_settings/20221027 enrichment_settings.csv")

file_structure <-
  create_folder_structure(dataset_id = dataset_metadata$dataset_id,
                          dataset_run_id = dataset_metadata$version,
                          main_folder = "results/Clustering_results",
                          folders = c("data" = "data", 
                                      "svg" = "svg",
                                      "UMAP" = "UMAP",
                                      "PCA" = "PCA",
                                      "treemap" = "svg/treemap",
                                      "heatmap" = "svg/heatmap",
                                      "distance" = "distance",
                                      "graph" = "graph",
                                      "enrichment" = "enrichment",
                                      "clustering" = "clustering"))
###

clustering_data <- 
  file_structure %>% 
  map(function(x) x$clustering) %>% 
  map(function(x) {
    file_ <- paste(x, "final_consensus.tsv", sep = "/")
    if(file.exists(file_)) {
      read_tsv(file_)
    } else {
      NULL
    }
  }) %>% 
  bind_rows(.id = "dataset_id") %>% 
  mutate(cluster = as.character(cluster))  %>% 
  filter(dataset_id %in% c("singlecell", "celline"))

pca_data <- 
  file_structure %>% 
  map(function(x) x$PCA) %>% 
  map(function(x) {
    file_ <- paste(x, "PCA.rds", sep = "/")
    if(file.exists(file_)) {
      readRDS(file_)
    } else {
      NULL
    }
  }) 


umap_data <- 
  file_structure %>% 
  map(function(x) x$UMAP) %>% 
  map(function(x) {
    file_ <- paste(x, "UMAP.tsv", sep = "/")
    if(file.exists(file_)) {
      read_tsv(file_)
    } else {
      NULL
    }
  }) %>% 
  bind_rows(.id = "dataset_id") %>% 
  filter(dataset_id %in% c("singlecell", "celline"))

umap_poly_data <- 
  file_structure %>% 
  map(function(x) x$UMAP) %>% 
  map(function(x) {
    file_ <- paste(x, "UMAP_polygons.tsv", sep = "/")
    if(file.exists(file_)) {
      read_tsv(file_)
    } else {
      NULL
    }
  }) %>% 
  bind_rows(.id = "dataset_id") %>% 
  mutate(cluster = as.character(cluster))  %>% 
  filter(dataset_id %in% c("singlecell", "celline"))

umap_centers <- 
  file_structure %>% 
  map(function(x) x$UMAP) %>% 
  map(function(x) {
    file_ <- paste(x, "cluster_centers.tsv", sep = "/")
    if(file.exists(file_)) {
      read_tsv(file_)
    } else {
      NULL
    }
  }) %>% 
  bind_rows(.id = "dataset_id") %>% 
  mutate(cluster = as.character(cluster))  %>% 
  filter(dataset_id %in% c("singlecell", "celline"))


cluster_annotation <-
  # list(tissue = read_csv("data/annotation/Tissue 20211013.csv"),
  #      blood = read_csv("data/annotation/Blood 20211013.csv"),
  #      celline = read_csv("data/annotation/Celline 20211013.csv"),
  #      brain = read_csv("data/annotation/Brain 20211013.csv")) %>% 
  list.files("data/annotation", full.names = T) %>% 
  {.[grepl("20221107", .)]} %>% 
  set_names(., gsub(".*\\d |\\.xlsx$", "", .)) %>% 
  map(. %>% read_excel()) %>% 
  bind_rows(.id = "dataset_id") %>%
  filter(dataset_id != "brain") %>% 
  rename(dataset_id = 1,
         cluster = 2,
         n_genes = 3,
         ann_specificity = 5, 
         ann_function = 6,
         reliability = 7,
         external_comment = 8,
         internal_comment = 9) %>% 
  mutate(cluster = as.character(cluster),
         full_annotation = paste(ann_specificity, 
                                 ann_function, 
                                 sep = " - "))

# ------ Nearest string function ----



complete_palette_neighbor <- 
  function(pal_terms, pal, force_colors = "") {
    
    pal_terms %>% 
      enframe("i", "name") %>% 
      select(-i) %>% 
      expand_grid(term = names(pal)) %>% 
      distinct() %>% 
      left_join(pal %>% 
                  enframe("term", "color"),
                by = c("term")) %>% 
      distinct() %>% 
      mutate(strsim = stringdist::stringsim(name, term)) %>% 
      arrange(-strsim) %>% 
      group_by(name) %>% 
      slice(1) %>% 
      ungroup() %>% 
      select(name, color) %>% 
      mutate(color = case_when(name %in% names(force_colors) ~ force_colors[match(name, names(force_colors))],
                               T ~ color)) %>% 
      deframe()
    
  }

# ------ Generate UMAP plot ----

umap_anno_centers <-
  umap_poly_data %>%
  filter(sub_cluster == 1) %>% 
  group_by(dataset_id, cluster) %>% 
  summarise(Xmax = max(X),
            Xmin = min(X),
            Ymax = max(Y),
            Ymin = min(Y),
            X = median(X),
            Y = median(Y)) %>% 
  ungroup() %>% 
  left_join(cluster_annotation %>% 
              select(dataset_id, 
                     cluster,
                     ann_specificity, 
                     ann_function) %>% 
              mutate(annotation = paste(ann_specificity, ann_function, 
                                        sep = "\n"))) 




# Label center adjustments: 
# cluster, x_nudge, y_nudge, hjust, vjust
label_adj <- 
  c(c(1, 0, 1, 0.5, 0),
    c(2, 1, 0, 0, 0),
    c(3, 1, 1, 0, 0),
    c(4, 0, 1, 0, 0),
    c(5, 1, 0, 0, 0.5),
    c(6, -1, 0, 1, 0),
    c(7, 0, -1, 0.5, 0.3),
    c(8, -1, 0, 1, 0.5),
    c(9, 0, -1, 0.5, 1),
    c(10, 0, -1, 0.5, 1),
    c(11, 0, -1, 0.5, 1),
    c(12, 1, 0, 0, 0.5),
    c(13, -1, 0, 1, 0.5),
    c(14, 1, 0, 0, 0.5),
    c(15, 1, 0, 0, 0),
    c(16, 1, 0, 0, 1),
    c(17, 0, -1, 0, 1),
    c(18, -1, 0, 1, 0.5),
    c(19, -1, 0, 1, 0.5),
    c(20, 0, 1, 0.5, 0),
    c(21, 0, -1, 0, 0.5),
    c(22, 0, 0, 0, 0),
    c(23, 1, 0, 0, 0),
    c(24, -1, 0, 1, 0.5),
    c(25, -1, 0, 1, 0.5),
    c(26, 0, -1, 0, 1),
    c(27, 1, 1, 0, 0),
    c(28, 1, 0, 0, 0.5),
    c(29, -1, 0, 1, 0.5),
    c(30, 0, -1, 0, 1),
    c(31, -1, 0, 1, 0.5),
    c(32, 0, -1, 0.5, 1),
    c(33, -1, 0, 1, 0.5),
    c(34, 1, 0, 0, 0.5),
    c(35, 0, 1, 0.5, 0),
    c(36, 1, 0, 0, 0),
    c(37, 0, 1, 0.5, 0),
    c(38, 0, -1, 0.5, 1),
    c(39, 0, 0, 1, 0),
    c(40, -1, 0, 1, 0),
    c(41, -1, 0, 1, -0.5),
    c(42, 0, -1, 0, 1),
    c(43, 1, 0, 0, 0.5),
    c(44, 0, -1, 0, 0),
    c(45, 1, 1, 0, 1),
    c(46, 1, 1, 1, 0),
    c(47, 1, 0, 0, 0.5),
    c(48, 0, -1, 0.5, 1),
    c(49, 0, 1, 0.5, 0),
    c(50, -1, 0, 1, 0.5),
    c(51, 1, 0, 0, 0.5),
    c(52, 0, -1, 0, -0.5),
    c(53, 1, 1, 0, 1),
    c(54, 0, -1, 0.5, 1),
    c(55, 0, -1, 0.5, 1),
    c(56, 1, 1, 0, 1),
    c(57, 0, 1, 0, 1),
    c(58, -1, 0, 1, 0.5),
    c(59, 0, -1, 0, 1),
    c(60, 1, 0, 0, 0.5),
    c(61, 1, 0, 0, 0.5),
    c(62, 1, 0, 0, 0.5),
    c(63, -1, 0, 1, 0.5),
    c(64, 0, -1, 0, 0),
    c(65, -1, 0, 1, 0),
    c(66, 0, 1, 1, 1),
    c(67, 0, 1, 1, 0),
    c(68, 0, 1, 1, 0)) %>% 
  matrix(ncol = 5, byrow = T) %>% 
  set_colnames(c("cluster", "x_nudge", "y_nudge", "hjust", "vjust")) %>% 
  as_tibble() %>% 
  mutate(cluster = as.character(cluster))

umap_anno_centers_adj <- 
  umap_anno_centers %>% 
  left_join(label_adj) %>% 
  group_by_all() %>% 
  mutate(X = c(Xmin, X, Xmax)[x_nudge + 2],
         Y = c(Ymin, Y, Ymax)[y_nudge + 2]) %>% 
  ungroup()



# ------ Generate specificity palettes ----

plot_pal <- 
  complete_palette_neighbor(pal_terms = unique(cluster_annotation$ann_specificity),
                            pal = region_palette,
                            force_colors = c("Non-specific" = "gray50"))



# ------ UMAP vanilla ----

umap_poly_data %>% 
  left_join(cluster_annotation %>% 
              select(cluster, ann_specificity)) %>% 
  ggplot(aes(X, Y)) +
  geom_point(data = umap_data,
             aes(UMAP_1_scaled,
                 UMAP_2_scaled),
             color = "darkgray",
             size = 0.1,
             inherit.aes = F) +
  geom_polygon(aes(group = polygon_id,
                   fill = as.numeric(cluster)),
               show.legend = F,
               color = NA,
               size = 0.1,
               alpha = 0.8) +
  facet_wrap(~dataset_id) +
  # scale_fill_identity() +
  scale_color_gradientn(colors = ggthemes::tableau_color_pal(palette = "Classic Cyclic")(13)) +
  scale_fill_gradientn(colors = ggthemes::tableau_color_pal(palette = "Classic Cyclic")(13)) +
  # scale_color_manual(values = plot_pal) +
  # scale_fill_manual(values = plot_pal) +
  coord_fixed() +
  theme_void() 
ggsave(savepath("UMAP basic.pdf"),
       width = 12, height = 8)
ggsave(savepath("UMAP basic.png"),
       width = 12, height = 8)


# ------ UMAP only numbers ----
umap_poly_data %>% 
  left_join(cluster_annotation %>% 
              select(cluster, ann_specificity)) %>% 
  ggplot(aes(X, Y)) +
  geom_point(data = umap_data,
             aes(UMAP_1_scaled,
                 UMAP_2_scaled),
             color = "darkgray",
             size = 0.1,
             inherit.aes = F) +
  geom_polygon(aes(group = polygon_id,
                   fill = ann_specificity),
               show.legend = F,
               color = NA,
               size = 0.1,
               alpha = 0.8) +
  geom_point(data = umap_anno_centers %>% 
               left_join(cluster_annotation %>% 
                           select(cluster, ann_specificity)),
             aes(X,
                 Y, 
                 color = ann_specificity),
             # color = "darkgray",
             show.legend = F,
             size = 5,
             inherit.aes = F) +
  geom_text(data = umap_anno_centers,
            aes(X,
                Y, 
                label = cluster),
            color = "white",
            size = 3,
            inherit.aes = F) +
  facet_wrap(~dataset_id) +
  # scale_fill_identity() +
  scale_color_manual(values = plot_pal) +
  scale_fill_manual(values = plot_pal) +
  coord_fixed() +
  theme_void() 
ggsave(savepath("UMAP only cluster number cellcolor.pdf"),
       width = 12, height = 8)
ggsave(savepath("UMAP only cluster number cellcolor.png"),
       width = 12, height = 8)

umap_poly_data %>% 
  left_join(cluster_annotation %>% 
              select(cluster, ann_specificity)) %>% 
  ggplot(aes(X, Y)) +
  geom_point(data = umap_data,
             aes(UMAP_1_scaled,
                 UMAP_2_scaled),
             color = "darkgray",
             size = 0.1,
             inherit.aes = F) +
  geom_polygon(aes(group = polygon_id,
                   fill = as.numeric(cluster)),
               show.legend = F,
               color = NA,
               size = 0.1,
               alpha = 0.8) +
  geom_point(data = umap_anno_centers %>% 
               left_join(cluster_annotation %>% 
                           select(cluster, ann_specificity)),
             aes(X,
                 Y, 
                 color = as.numeric(cluster)),
             # color = "darkgray",
             show.legend = F,
             size = 5,
             inherit.aes = F) +
  geom_text(data = umap_anno_centers,
            aes(X,
                Y, 
                label = as.numeric(cluster)),
            color = "white",
            size = 3,
            inherit.aes = F) +
  facet_wrap(~dataset_id) +
  # scale_fill_identity() +
  scale_color_gradientn(colors = ggthemes::tableau_color_pal(palette = "Classic Cyclic")(13)) +
  scale_fill_gradientn(colors = ggthemes::tableau_color_pal(palette = "Classic Cyclic")(13)) +
  # scale_color_manual(values = plot_pal) +
  # scale_fill_manual(values = plot_pal) +
  coord_fixed() +
  theme_void() 
ggsave(savepath("UMAP only cluster number rancolor.pdf"),
       width = 12, height = 8)
ggsave(savepath("UMAP only cluster number rancolor.png"),
       width = 12, height = 8)



# ------ UMAP annotation ----
# Generate palette from specificity annotation


umap_plot_original <- 
  umap_poly_data %>% 
  left_join(cluster_annotation %>% 
              select(dataset_id, cluster, ann_specificity)) %>% 
  ggplot(aes(X, Y)) +
  geom_point(data = umap_data,
             aes(UMAP_1_scaled,
                 UMAP_2_scaled),
             color = "darkgray",
             size = 0.1,
             inherit.aes = F) +
  geom_polygon(aes(group = polygon_id,
                   fill = ann_specificity),
               show.legend = F,
               color = NA,
               size = 0.1,
               alpha = 0.8) +
  geom_point(data = umap_anno_centers %>% 
               left_join(cluster_annotation %>% 
                           select(dataset_id, cluster, ann_specificity)),
             aes(X,
                 Y, 
                 color = ann_specificity),
             # color = "darkgray",
             show.legend = F,
             size = 5,
             inherit.aes = F) +
  geom_text(data = umap_anno_centers,
            aes(X,
                Y, 
                label = as.numeric(cluster)),
            color = "white",
            size = 3,
            inherit.aes = F) +
  
  geom_text(data = umap_anno_centers,
            aes(label = annotation),
            lineheight = 0.8, 
            size = 2) +
  facet_wrap(~dataset_id) +
  # scale_fill_identity() +
  scale_color_manual(values = plot_pal) +
  scale_fill_manual(values = plot_pal) +
  coord_fixed() +
  theme_void() 
ggsave(savepath("UMAP annotation original.pdf"),
       plot = umap_plot_original,
       width = 16, height = 8)
ggsave(savepath("UMAP annotation original.png"),
       plot = umap_plot_original,
       width = 16, height = 8)


# ------ annotation stats ----


anno_joined$reliability_workshop %>% unique
anno_joined %>%
  filter(!is.na(specificity_workshop)) %>% 
  filter(!is.na(reliability_workshop)) %>% 
  mutate(reliability_original = factor(str_to_sentence(reliability_original),
                                       c("Very low", "Low", "Medium", "High")),
         reliability_workshop = reliability_workshop %>% 
           str_remove("\\?") %>% 
           str_to_sentence() %>% 
           factor(c("Very low", "Low", "Medium", "High"))) %>%
  filter(!is.na(reliability_workshop)) %>% 
  group_by(reliability_original, reliability_workshop) %>% 
  count() %>% 
  ggplot(aes(reliability_original, reliability_workshop, 
             fill = n, label = n)) +
  geom_tile() +
  geom_text() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, name = "YlOrRd")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_fixed() +
  labs(title = "Reliability",
       x = "HPA",
       y = "Workshop")

ggsave("results/reliability matrix.pdf",
       width = 4, height = 4)



anno_joined %>% 
  select(cluster, HPA = specificity_original, Workshop = specificity_workshop) %>% 
  gather(Type, Annotation, HPA, Workshop) %>% 
  mutate(Specificity = case_when(tolower(Annotation) == "non-specific" ~ "Non-specific",
                                 is.na(Annotation) ~ "Not annotated",
                                 T ~ "Specific")) %>% 
  group_by(Type, Specificity) %>% 
  count() %>% 
  ggplot(aes(n, Type, fill = Specificity)) +
  geom_col() +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5),
            color = "white") +
  scale_fill_manual(values = c("Non-specific" = "gray",
                               "Not annotated" = "gray30",
                               "Specific" = "orangered")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid = element_blank())
ggsave("results/number annotated as specific.pdf",
       width = 5, height = 2)



anno_joined %>% 
  filter(!is.na(specificity_workshop)) %>% 
  select(cluster, HPA = specificity_original, Workshop = specificity_workshop) %>% 
  gather(Type, Annotation, HPA, Workshop) %>% 
  mutate(Specificity = case_when(tolower(Annotation) == "non-specific" ~ "Non-specific",
                                 is.na(Annotation) ~ "Not annotated",
                                 T ~ "Specific")) %>% 
  group_by(Type, Specificity) %>% 
  count() %>% 
  ggplot(aes(n, Type, fill = Specificity)) +
  geom_col() +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5),
            color = "white") +
  scale_fill_manual(values = c("Non-specific" = "gray",
                               "Not annotated" = "gray30",
                               "Specific" = "orangered")) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid = element_blank())
ggsave("results/number annotated as specific only annotated.pdf",
       width = 5, height = 2)

