
library(tidyverse)
library(ggtern)
source("scripts/functions_utility.R")


plot_data <- 
  read_tsv("data/processed/sc exp sum.tsv")

consensus_colors <- 
  read_tsv("data/colors/colors_consensus.tsv")

consensus_palette <- 
  consensus_colors %>% 
  select(sample, color) %>% 
  deframe()

plot_data %>%
  select(-ntpm) %>% 
  spread(class, frac) %>% 
  ggtern(aes(secreted, membrane, intracellular,
             color = cell_type)) +
  geom_point(show.legend = F) +
  geom_text(data = . %>% 
              filter(secreted > 0.2 |
                       membrane > 0.5 |
                       intracellular > 0.86),
            aes(label = cell_type),
            show.legend = F, 
            size = 3,
            vjust = 1,
            hjust = 1.05) +
  theme_bw() +
  # xlab("More") +                       #replace default axis labels
  # ylab("More membrane") +
  # zlab("Blue") +
  scale_color_manual(values = consensus_palette)

ggsave(savepath("protein loc frac tern.pdf"),
       width = 6, height = 6)
