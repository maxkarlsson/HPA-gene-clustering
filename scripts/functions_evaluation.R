
# Evaluation 

eval <- function(clustering) {
  con_index <-
    connectivity(clustering)
  dunn <- 
    dunn()
  
  bio_index <- 
    if(require("Biobase") && require("annotate") && require("GO.db") &&
     require("moe430a.db")) {
    BHI(cluster, annotation="moe430a.db", names=rownames(clustering), category="all")
  }
  
  list(connectivity = con_index,
       dunn = dunn_index,
       BHI = bio_index)
}