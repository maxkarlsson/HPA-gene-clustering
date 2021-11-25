
library(tidyverse)

raw_data <-
  readRDS("data/processed/combined_HPA_expression_data.RDS")

color_celline <- 
  read_tsv("data/colors/Cell lines.txt",
           col_names = c("organ", "sample", "color"))
color_singlecell <- 
  read_tsv("data/colors/Single cell type.txt",
           col_names = c("organ", "sample", "color"))
color_tissue <- 
  read_tsv("data/colors/Tissues.txt",
           col_names = c("organ", "sample", "color"))

all_colors <- 
  bind_rows(celline = color_celline,
            singlecell = color_singlecell,
            tissue = color_tissue,
            .id = "dataset_id") %>% 
  bind_rows(.,
            filter(.,
                   sample %in% c("endometrium 1",
                                 "stomach 1",
                                 "skin 1",
                                 "lymphoid tissue")) %>% 
              mutate(sample = case_when(sample == "endometrium 1" ~ "endometrium",
                                        sample == "stomach 1" ~ "stomach",
                                        sample == "skin 1" ~ "skin",
                                        sample == "lymphoid tissue" ~ "lymphoid system")))


samples_in_data <-
  raw_data %>% 
  map(. %>% 
        select(-1) %>% 
        names() %>% 
        enframe("i", "sample")) %>% 
  bind_rows(.id = "dataset_id") %>% 
  select(-i) %>% 
  mutate(sample = ifelse(dataset_id == "singlecell_sample", 
                         str_extract(sample, "_.*_") %>% 
                           gsub("_", "", .),
                         gsub("_.*", "", sample))) %>% 
  distinct()



### Consensus colors

consensus_colors <- 
  samples_in_data %>% 
  filter(grepl("consensus", dataset_id)) %>% 
  mutate(dataset_id = gsub("_consensus", "", dataset_id)) %>% 
  left_join(all_colors %>% 
              select(sample, color) %>% 
              distinct(),
            by = "sample") %>%
  group_by(sample) %>% 
  mutate(n = n_distinct(color)) %>% 
  ungroup() %>% 
  # filter(n > 1) %>% 
  filter(!(n > 1 & 
             ((dataset_id == "singlecell" &
                 color != "#b30000") |
                (sample == "choroid plexus" &
                color != "#A7DACD") |
                (sample == "pituitary gland" &
                   color != "#7F6A9C") |
                (sample == "retina" &
                   color != "#FFEF78")))) %>% 
  # mutate(case) %>% 
  select(dataset_id, sample, color)


consensus_colors %>% 
  arrange(color) %>% 
  group_by(dataset_id) %>% 
  mutate(y = row_number(),
         x = y %% 8) %>% 
  ggplot(aes(x, y, label = sample, fill = color)) +
  geom_point(shape = 22,
             size = 7) +
  geom_text(size = 2) +
  facet_wrap(~dataset_id) +
  scale_fill_identity() +
  theme_bw()

consensus_colors %>% 
  write_tsv("data/colors/colors_consensus.tsv")


## region / sample colors

region_colors <- 
  samples_in_data %>% 
  filter(!(grepl("consensus", dataset_id) &
             dataset_id != "celline_consensus")) %>% 
  filter(!grepl("brain", dataset_id)) %>% 
  left_join(all_colors %>% 
              select(sample, organ, color) %>% 
              distinct(),
            by = "sample") %>%
  group_by(sample) %>% 
  mutate(n = n_distinct(color)) %>% 
  ungroup() %>% 
  # filter(n > 1) %>% 
  # filter(sample == "choroid plexus") %>% 
  # filter(!(sample == "choroid plexus" &
  #            color != "#A7DACD"))
  group_by(sample) %>% 

  filter(!(n > 1 & 
             ((any(organ == "Brain") &
                 dataset_id == "tissue_region" & 
                 organ != "Brain") |
                (dataset_id == "singlecell_sample" &
                   color != "#b30000") |
                (sample == "choroid plexus" &
                   color != "#A7DACD") |
                (sample == "pituitary gland" &
                   color != "#7F6A9C") |
                (sample == "retina" &
                   color != "#FFEF78")))) %>% 
  # mutate(case) %>% 
  ungroup() %>% 
  select(dataset_id, sample, color)


region_colors %>% 
  arrange(color) %>% 
  group_by(dataset_id) %>% 
  mutate(y = row_number(),
         x = y %% 8) %>% 
  ggplot(aes(x, y, label = sample, fill = color)) +
  geom_point(shape = 22,
             size = 7) +
  geom_text(size = 2) +
  facet_wrap(~dataset_id) +
  scale_fill_identity() +
  theme_bw()

region_colors %>% 
  write_tsv("data/colors/colors_regions.tsv")



