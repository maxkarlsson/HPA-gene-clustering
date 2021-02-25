
# Create random clustering

random_cluster <- function (k,l, names) {
  cluster <- sample(1:k, l, replace = T)
  return(data.frame(gene = names, cluster = cluster))
}


# Evaluation 

eval <- function(clustering, genes, dist, id = NULL) { # add distance
  
  #clustering <- dplyr::rename(clustering, cluster = value)
  
  r <- 
    clustering %>%
    left_join(genes, by = c("gene" = "enssscg_id")) %>%
    na.omit() %>%
    dplyr::select(-gene) 
  
  
  clusters <- as.numeric(as.character(r$cluster))
  
  con_index <-
    connectivity(distance = dist, cluster = as.numeric(as.character(clustering$cluster)))
  
  dunn_index <- 
    dunn(distance = dist, cluster = as.numeric(as.character(clustering$cluster)))
  
  bio_index <- 
    if(require("Biobase") && require("annotate") && require("GO.db") &&
       require("org.Hs.eg.db")) {
      BHI(clusters, annotation="org.Hs.eg.db", names=r$entrez, category="all") 
    }

  
  if(!is.null(id)) {
    return(list(connectivity = con_index,
                dunn = dunn_index,
                BHI = bio_index,
                id = id))
  }
  else {
    return(list(connectivity = con_index,
                dunn = dunn_index,
                BHI = bio_index))
  }

}