library(tidyverse)
library(NOISeq)
library(uwot)
library(pcaMethods)
library(magrittr)
library(factoextra)
library(cluster)
library(kohonen)
library(pheatmap)
library(clValid)
library(ggalt)
library(sf)
library(concaveman)
library(ggthemes)
library(ggrepel)
#library(multidplyr) #Patition a data frame across multiple worker processes in to provide simple multicore parallelism.
#library(FNN) #Fast Nearest Neighbor Search Algorithms and Applications
#library(NbClust) #provides 30 indices for determining the number of clusters and proposes to user the best clustering scheme from the different results obtained by varying all combinations of number of clusters, distance measures, and clustering methods.

data <-
  read_delim("Data/pig_filtered_sample_data.tab", delim = "\t")
```