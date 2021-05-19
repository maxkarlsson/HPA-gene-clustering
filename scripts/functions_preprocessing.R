
# pTPM normalization function

ptpm_normalize <- function(df, sample_col, tpm_col, format = "wide") {
  df_normalized <-
    df %>% 
    group_by(.data[[sample_col]]) %>%
    mutate(sum_tpm = sum(.data[[tpm_col]]),
           ptpm = .data[[tpm_col]]* 1e6 / sum_tpm) %>%
    select(-sum_tpm) %>%
    ungroup() 
  
  if (format == "wide") {
    df_normalized <-
      df_normalized %>%
      select(1, 2, ptpm) %>%
      spread(key = .data[[sample_col]], value = ptpm) %>%
      column_to_rownames("enssscg_id")
    return(df_normalized)
  }
  
  else if (format == "long") {
    return(df_normalized)
  }
}


# TMM normalization function

tmm_normalize <- function(x, ...) {
  median_column <-
    apply(x,
          MARGIN = 1,
          median)
  x[, dim(x)[2] + 1] <- median_column
  norm_data <- NOISeq::tmm(x, refColumn = dim(x)[2], ...)
  norm_data[, -(dim(x)[2])] 
}


# Data scaling function

data_scaling <- function(df,col_value,col_gene,col_sample, 
                         log_transform = F, m = "zscore") {
  
  if (log_transform == T) {
    df <-
      df %>%
      mutate(tmm = log(.data[[col_value]] + 1))
  } else {
    df <-
      df %>%
      mutate(tmm = .data[[col_value]])
  }
  
  if (m == "zscore") { 
    df_scaled <-
      df %>%
      group_by(.data[[col_gene]]) %>%
      mutate(mean = mean(tmm),
             sd = sd(tmm)) %>%
      ungroup() %>%
      mutate(z_score = (tmm - mean) / sd) %>%
      select(.data[[col_gene]], .data[[col_sample]], z_score) %>% 
      mutate(z_score = ifelse(is.na(z_score), 0, z_score)) %>%
      spread(.data[[col_sample]], z_score) %>% 
      column_to_rownames(col_gene)
  }
  
  else if (m == "min-max") {
    df_scaled <-
      df %>% 
      group_by(.data[[col_gene]]) %>%
      mutate(min = min(tmm),
             max = max(tmm)) %>%
      ungroup() %>%
      mutate (min_max = (tmm-min) / (max-min)) %>%
      select(.data[[col_gene]], .data[[col_sample]], min_max) %>% 
      mutate(min_max = ifelse(is.na(min_max), 0, min_max)) %>%
      spread(.data[[col_sample]], min_max) %>% 
      column_to_rownames(col_gene)
  }
  
  
  else if (m == "max") {
    df_scaled <-
      df %>% 
      group_by(.data[[col_gene]]) %>%
      mutate(max = max(tmm)) %>%
      ungroup() %>%
      mutate (max_scaled = tmm/max) %>%
      select(.data[[col_gene]], .data[[col_sample]], max_scaled) %>% 
      mutate(max_scaled = ifelse(is.na(max_scaled), 0, max_scaled)) %>%
      spread(.data[[col_sample]], max_scaled) %>% 
      column_to_rownames(col_gene)
    
  }
  
  list(method = m, scaled = df_scaled %>% as_tibble(rownames= col_gene))
}
