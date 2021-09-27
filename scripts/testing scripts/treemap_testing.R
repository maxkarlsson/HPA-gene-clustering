


library(tidyverse)
library(magrittr)
library(rrvgo)
library(treemapify)
library(pbapply)

cp_res_all <- 
  read_csv("results/Annotation reports/fischer test results/enrichment_results.csv")

ontology_pal <- 
  c("MF" = "#FF6F00",
    "CC" = "#C71B00",
    "BP" = "#018EA0")

GO_sim <-
  pblapply(c("BP" = "GO_BP",
             "CC" = "GO_CC", 
             "MF" = "GO_MF"),
           function(ont_) {
             GO_terms <-
               cp_res_all %>% 
               filter(test_type == ont_) %>% 
               pull(ID) %>% 
               unique()
             
             calculateSimMatrix(GO_terms, 
                                orgdb = "org.Hs.eg.db",
                                ont = str_extract(ont_, "BP|CC|MF"),
                                method = "Wang")
           })



plot_legend <-
  enframe(ontology_pal,
          "Ontology", 
          "Color") %>% 
  ggplot(aes(1, 1, fill = Ontology)) +
  geom_tile() +
  scale_fill_manual(values = ontology_pal)

ggpubr::get_legend(plot_legend) %>% 
  ggpubr::as_ggplot()
ggsave("results/Ontology legend.pdf",
       width = 1, height = 1)


plot_data <- 
  cp_res_all %>%
  filter(test_type %in% c("GO_BP",
                          "GO_CC", 
                          "GO_MF")) %>% 
  mutate(GeneFrac = as.numeric(gsub("\\/.*", "", GeneRatio)) /
           as.numeric(gsub(".*\\/", "", GeneRatio)),
         ontology = str_extract(test_type, "BP|CC|MF")) %>% 
  group_by(ontology, test_type, dataset_id, cluster) %>% 
  filter(n_distinct(ID) >= 2) %>% 
  group_by(ontology, test_type, dataset_id, cluster) %>% 
  do({
    g_data <<- .
    
    sim_mat <- 
      GO_sim[[unique(g_data$ontology)]][g_data$ID,
                                        g_data$ID]
    
    sim_mat %>%
      reduceSimMatrix(setNames(-log10(g_data$p.adjust), g_data$ID),
                      threshold=0.7,
                      orgdb="org.Hs.eg.db") %>%
      as_tibble() %>%
      rename(term_cluster = cluster)
    
    
  }) %>% 
  
  ungroup()


spread_colors <- 
  function(color, n, colorrange = 0.5) {
    
    if(n == 1) {
      return(color)
    } else {
      colorRamp(c("white", color, "black"))(seq(0.5 - colorrange / 2,
                                                0.5 + colorrange / 2, 
                                                length.out = n)) %>% 
        as_tibble() %$% 
        rgb(V1, V2, V3, maxColorValue = 255) 
      
    }
    
  }


plot_settings <-
  plot_data %>% 
  select(dataset_id, cluster) %>% 
  distinct() %>% 
  mutate(i = row_number())

plots <- 
  pblapply(plot_settings$i,
           function(i) {
             a %>% 
               filter(dataset_id == plot_settings$dataset_id[i],
                      cluster == plot_settings$cluster[i]) %>% 
               left_join(enframe(ontology_pal,
                                 "ontology", 
                                 "color"),
                         by = "ontology") %>% 
               group_by(parent) %>% 
               mutate(nn = n_distinct(go),
                      color2 = spread_colors(unique(color), n_distinct(go), colorrange = 0.2)) %>% 
               ggplot(aes(area = score, subgroup = parentTerm)) +
               
               geom_treemap(aes(fill = color2), 
                            color = "black",
                            show.legend = T) +
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
               
               scale_fill_identity() +
               theme_void() +
               ggtitle(paste(plot_settings$dataset_id[i],
                             "cluster",
                             plot_settings$cluster[i]))
           })


pdf("results/GO treemaps.pdf")
plots
dev.off()

