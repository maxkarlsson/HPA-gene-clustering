

# Get gene class

gene_classification <- 
  enrichment_settings %>% 
  filter(grepl("specificity", id)) %>% 
  group_by(id) %>% 
  do({
    g_data <<- .
    
    get_specificity_db(enrichment_settings, g_data$id)
    
  }) %>% 
  ungroup()

# Get long gene class

## Alt 1

specificity_long <- 
  enrichment_settings %>% 
  filter(grepl("specificity", id)) %>% 
  group_by(id) %>% 
  do({
    g_data <<- .
    
    get_db(enrichment_settings, g_data$id)
    
    
  }) %>% 
  ungroup()

## Alt 2
specificity_long2 <- 
  gene_classification %>% 
  mutate(enhanced_tissues = gsub(", ", ";", enhanced_tissues)) %>% 
  separate_rows(enhanced_tissues, sep = ",") %>% 
  mutate(enhanced_tissues = gsub(";", ", ", enhanced_tissues)) 
  

