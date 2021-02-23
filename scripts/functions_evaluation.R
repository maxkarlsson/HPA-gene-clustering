
# Create random clustering


# Evaluation 

eval <- function(clustering, genes, distance, id = NULL) { # add distance
  dplyr::rename(clustering, cluster = value)
  
  r <- 
    clustering %>%
    left_join(genes, by = c("gene" = "enssscg_id")) %>%
    na.omit() %>%
    dplyr::select(-gene) %>%
  
  
  clusters <- as.numeric(as.character(r$cluster))
  
  con_index <-
    connectivity(distance = distance, cluster = clustering$cluster)
  
  dunn <- 
    dunn(distance = distance, cluster = clustering$cluster)
  
  bio_index <- 
    if(require("Biobase") && require("annotate") && require("GO.db") &&
       require("org.Hs.eg.db")) {
      BHI(clusters, annotation="org.Hs.eg.db", names=r$entrez, category="all") 
    }
  
  list(connectivity = con_index,
       dunn = dunn_index,
       BHI = bio_index)
}