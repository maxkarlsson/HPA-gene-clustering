




library(tidyverse)




res_clusters <- read_csv("data/processed/num_clusters.csv")


cluster_range <- 
  tibble(k = c(30, 110))


k_res_models <- 
  res_clusters %>% 
  pull(dataset_id) %>% 
  unique() %>% 
  set_names(., .) %>% 
  lapply(function(id_) lm(resolution ~ k, 
                          data = res_clusters %>% 
                            filter(dataset_id == id_)))

model_coefficients <-
  k_res_models %>% 
  map(coefficients) %>% 
  map(enframe) %>% 
  bind_rows(.id = "dataset_id") %>% 
  mutate(name = gsub("\\(|\\)", "", name)) %>% 
  spread(name, value)

predicted_resolutions <- 
  k_res_models %>% 
  map(function(model) predict(model, newdata = cluster_range) %>% 
        enframe("type", "resolution") %>% 
        mutate(type = ifelse(type == 1, 
                             "min", "max"))) %>% 
  bind_rows(.id = "dataset_id")



res_clusters %>% 
  ggplot(aes(k, resolution, color = dataset_id)) +
  geom_point() +
  # geom_smooth(method = "lm") +
  geom_abline(data = model_coefficients,
              aes(slope = k, intercept = Intercept)) +
  geom_hline(data = predicted_resolutions,
             aes(yintercept = resolution),
             linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = cluster_range$k, 
             linetype = "dashed") +
  facet_wrap(~dataset_id) +
  theme_bw()

range(predicted_resolutions$resolution)

predicted_resolutions %>% 
  spread(type, resolution)


###########


k_res_models_inverse <- 
  res_clusters %>% 
  pull(dataset_id) %>% 
  unique() %>% 
  set_names(., .) %>% 
  lapply(function(id_) lm(k ~ resolution, 
                          data = res_clusters %>% 
                            filter(dataset_id == id_)))

model_coefficients_inverse <-
  k_res_models_inverse %>% 
  map(coefficients) %>% 
  map(enframe) %>% 
  bind_rows(.id = "dataset_id") %>% 
  mutate(name = gsub("\\(|\\)", "", name)) %>% 
  spread(name, value)


model_coefficients_inverse$resolution %>% 
  max() %>% 
  {1/.}

seq(0.1, 13, by = 0.1) %>% 
  length()


