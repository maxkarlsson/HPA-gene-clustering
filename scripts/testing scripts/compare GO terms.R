


GO_terms <- 
  enrichment_settings %>% 
  filter(grepl("^GO", id)) %>% 
  filter(!is.na(file)) %>% 
  group_by(id) %>% 
  do({
    g_data <<- .
    

    
    raw_db <- 
      read_tsv(paste0("annotation_databases/", g_data$file))
    
    
    raw_db %>% 
      select(ensg_id = `Gene stable ID`, term = `GO term name`,     term_id = `GO term accession`) %>% 
      filter(!is.na(term)) %>% 
      distinct()
    
  }) %>% 
  ungroup()

#############
GO_terms <-
  database_terms %>%
  filter(grepl("^GO", id))


IT_terms <- read_csv("data/temp/HPAv21-gene-GO-211110.csv")

trans_terms <- 
  IT_terms %>% 
  group_by(ensembl_gene_id, ensembl_transcript_id) %>% 
  summarise(term = paste(sort(unique(external_id)), collapse = ";")) %>% 
  ungroup()

trans_terms

IT_gene_terms <- 
  IT_terms %>% 
  select(1, 3) %>% 
  distinct() %>% 
  filter(ensembl_gene_id %in% clustering_data$gene) %>% 
  group_by(external_id) %>% 
  filter(n_distinct(ensembl_gene_id) >= 8) %>% 
  ungroup()


n_distinct(IT_gene_terms$ensembl_gene_id)
n_distinct(GO_terms$ensg_id)

n_distinct(clustering_data$gene)

IT_gene_terms %>% 
  select(external_id) %>% 
  distinct() %>% 
  filter(!external_id %in% GO_terms$term_id)

GO_terms %>% 
  select(term_id) %>% 
  distinct() %>% 
  filter(!term_id %in% IT_gene_terms$external_id)
