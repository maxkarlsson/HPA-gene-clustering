


savepath <- 
  function(savename) { 
    result_folder <- paste0("results/", Sys.Date())
    dir.create(result_folder, showWarnings = FALSE)
    
    savename <-
      paste0(result_folder, "/", savename)
    
    
    return(savename)
    
  }

mix_colors <- 
  function(color1, color2, mix = 0) {
    sapply(1:length(color1),
           function(i) {
             colorRamp(c(color1[i], color2[i]))(mix[i]) %>% 
               rgb(maxColorValue = 255)
           }) 
  }

saveWidgetFix <- function (widget,file,...) {
  ## A wrapper to saveWidget which compensates for arguable BUG in
  ## saveWidget which requires `file` to be in current working
  ## directory.
  wd<-getwd()
  on.exit(setwd(wd))
  outDir<-dirname(file)
  file<-basename(file)
  setwd(outDir);
  saveWidget(widget,file=file,...)
}

aggregate_data <- 
  function(data, meta_data_, consensus_hierarchy, saveid, process_function) {
    cat(paste0(saveid, "\n"))
    savefile <- paste0("data/processed/", saveid, ".rds")
    if(file.exists(savefile)) {
      cat("File already present.\n")
      outdata <- 
        read_rds(savefile)
    } else {
      cat("Generating files.\n")
      
      data <- 
        data %>% 
        process_function()
      
      tissue_dataset_data <- 
        meta_data_ %>%  
        # filter(scilifelab_id %in% run_samples$scilifelab_id) %>% 
        
        group_by(dataset, platform, region) %>% 
        do({
          sample_ids <<- 
            .$scilifelab_id
          
          data %>%
            select(ensg_id, all_of(sample_ids)) %>%
            column_to_rownames("ensg_id") %>% 
            apply(MARGIN = 1,
                  function(x) mean(x, na.rm = T)) %>% 
            enframe("ensg_id", "tmm")
        }) %>%
        ungroup()
      
      tissue_dataset_sd <-
        meta_data_ %>%  
        # filter(norm_batch == "tissue") %>% 
        # filter(scilifelab_id %in% run_samples$scilifelab_id) %>% 
        
        group_by(dataset, platform, region) %>% 
        do({
          sample_ids <<- 
            .$scilifelab_id
          
          data %>%
            select(ensg_id, all_of(sample_ids)) %>%
            column_to_rownames("ensg_id") %>% 
            apply(MARGIN = 1,
                  function(x) sd(x, na.rm = T)) %>% 
            enframe("ensg_id", "sd")
        }) %>%
        ungroup()
      
      tissue_dataset_sd_log <-
        meta_data_ %>%  
        # filter(norm_batch == "tissue") %>% 
        # filter(scilifelab_id %in% run_samples$scilifelab_id) %>% 
        
        group_by(dataset, platform, region) %>% 
        do({
          sample_ids <<- 
            .$scilifelab_id
          
          data %>%
            select(ensg_id, all_of(sample_ids)) %>%
            column_to_rownames("ensg_id") %>% 
            apply(MARGIN = 1,
                  function(x) sd(log10(x + 1), na.rm = T)) %>% 
            enframe("ensg_id", "sd")
        }) %>%
        ungroup()
      
      tissue_data <- 
        tissue_dataset_data %>% 
        group_by(region, ensg_id) %>% 
        top_n(1, tmm) %>% 
        slice(1) %>% 
        ungroup()
      
      consensus_data <-
        tissue_data %>% 
        mutate(region = tolower(region)) %>%
        left_join(consensus_hierarchy %>% 
                    select(content_name, consensus = consensus_content_name),
                  by = c("region" = "content_name")) %>%
        group_by(consensus, ensg_id) %>%
        top_n(1, tmm) %>% 
        slice(1) %>% 
        ungroup()
      
      
      outdata <- 
        list(tissue_dataset_data = tissue_dataset_data, 
             tissue_data = tissue_data, 
             consensus_data = consensus_data,
             tissue_dataset_sd = tissue_dataset_sd,
             tissue_dataset_sd_log = tissue_dataset_sd_log)
      
      saveRDS(outdata, savefile)
      
      
    }
    cat("Done.\n")
    return(outdata)
  }



compile_example_data <- 
  function(file_name) {
    
    brain_samples <- 
      meta_data %>% 
      filter(dataset %in% c("hpa brain", "pfc brain"))
    
    
    all_data <- 
      list(tissue = 
             readRDS("../HPA-normalization2021/data/processed/aggregated_final_tissue_data2.rds")$tissue_data %>% 
             select(ensg_id, tmm, region) %>% 
             spread(region, tmm),
           brain = 
             vroom("../../Data/HPA/normalized/tissue_tmmnorm.tsv") %>% 
             select(1, all_of(brain_samples$scilifelab_id)) %>% 
             gather(scilifelab_id, tmm, -1) %>% 
             left_join(brain_samples) %>% 
             group_by(sample, ensg_id) %>% 
             summarise(tmm = mean(tmm)) %>% 
             select(ensg_id, tmm, sample) %>% 
             spread(sample, tmm),
           blood = 
             vroom("../../Data/HPA/normalized/blood_tmmnorm.tsv"),
           blood_consensus = 
             vroom("../../Data/HPA/normalized/blood_tmmnorm.tsv") %>% 
             gather(scilifelab_id, tmm, -1) %>% 
             left_join(meta_data) %>% 
             group_by(sample, ensg_id) %>% 
             summarise(tmm = mean(tmm)) %>% 
             select(ensg_id, tmm, sample) %>% 
             spread(sample, tmm),
           celline = 
             vroom("../../Data/HPA/normalized/celline_tmmnorm.tsv"),
           celline_consensus = 
             vroom("../../Data/HPA/normalized/celline_tmmnorm.tsv") %>% 
             gather(scilifelab_id, tmm, -1) %>% 
             left_join(meta_data) %>% 
             group_by(sample, ensg_id) %>% 
             summarise(tmm = mean(tmm)) %>% 
             select(ensg_id, tmm, sample) %>% 
             spread(sample, tmm),
           singlecell = 
             {
               single_cell_data <- 
                 vroom("data/meta/rna_single_cell_type_tissue.tsv") %>%
                 unite(id, Tissue, Cluster) %>% 
                 select(id, ensg_id = Gene, ptpm = pTPM) %>% 
                 spread(id, ptpm) %>% 
                 column_to_rownames("ensg_id") %>% 
                 as.matrix() 
               
               sc_tmm_factors <- 
                 
                 single_cell_data %>% 
                 
                 calc_Normfactors(method = "TMM", 
                                  
                                  # Use median column as reference distribution:
                                  refColumn = "median",
                                  
                                  # Trim parameters:
                                  logratioTrim = 0.3, 
                                  sumTrim = 0.3, 
                                  
                                  # Weighting should be done only if count data:
                                  doWeighting = F)
               
               
               
               single_cell_data <- t(t(single_cell_data) / sc_tmm_factors)
               
               single_cell_data %>% 
                 as_tibble(rownames = "ensg_id")
             }
      )
    
    saveRDS(all_data, paste0("data/processed/", file_name))
    
    return(all_data)
    
  }


calc_Normfactors <- 
  function (object, method = c("TMM", "quantile"), refColumn = NULL, 
            logratioTrim = 0.3, sumTrim = 0.05, doWeighting = TRUE, 
            Acutoff = -1e+10, quantile = 0.75) {
    method <- match.arg(method)
    if (is.matrix(object)) {
      if (is.null(refColumn)) 
        refColumn <- 1
      data <- object
      libsize <- colSums(data)
    } else {
      stop("calcNormFactors() only operates on 'matrix' objects")
    }
    
    if(refColumn == "median") {
      ref <- 
        apply(data, MARGIN = 1, median)
    } else {
      ref <- data[, refColumn]
    }
    
    f <- switch(method, TMM = apply(data, 2, NOISeq:::.calcFactorWeighted, 
                                    ref = ref, logratioTrim = logratioTrim, 
                                    sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff), 
                quantile = NOISeq:::.calcFactorQuantile(data, libsize, q = quantile))
    f <- f/exp(mean(log(f)))
    return(f)
  }



gene_fischer_test <- 
  function(cluster_membership, term_membership, alternative = "greater") {
    
    n_genes <- 
      term_membership %>% 
      unlist() %>% 
      n_distinct()
    
    fisher_data <-
      cluster_membership %>% 
      enframe("gene", "cluster") %>% 
      full_join(term_membership %>%
                  map(enframe, "row", "gene") %>% 
                  bind_rows(.id = "term")) %>% 
      mutate(term = factor(term),
             cluster = factor(cluster)) %>%
      group_by(cluster, term, .drop = F) %>% 
      count() %>%
      group_by(term) %>% 
      mutate(total_term = sum(n),
             total_nonterm = n_genes - total_term) %>% 
      filter(!is.na(cluster) &
               !is.na(term)) %>%
      group_by(cluster) %>% 
      mutate(total_in_cluster = sum(n),
             total_not_in_cluster = n_genes - total_in_cluster ) %>% 
      ungroup() 
    
    
    
    fisher_data %>% 
      group_by(cluster, term) %>% 
      do({
        g_data <<- .
        
        contingency_matrix <-
          matrix(c(g_data$n, 
                   g_data$total_in_cluster - g_data$n,
                   g_data$total_term - g_data$n,  
                   g_data$total_nonterm - (g_data$total_term - g_data$n)), 
                 nrow = 2,
                 dimnames =
                   list(c("Term", "Nonterm"),
                        c("In cluster", "Not in cluster")))
        
        contingency_matrix %>% 
          fisher.test(alternative = alternative) %>% 
          tidy() %>%
          mutate(term_overlap = paste(contingency_matrix[1], 
                                      contingency_matrix[1] + contingency_matrix[3],
                                      sep = "/"),
                 cluster_membership = paste(contingency_matrix[1], 
                                            contingency_matrix[1] + contingency_matrix[2],
                                            sep = "/"))
      }) %>% 
      ungroup() %>%
      mutate(p.value = ifelse(p.value == 0, .Machine$double.xmin, p.value),
             adj_p = p.adjust(p.value, method = "BH")) %>% 
      rename(odds_ratio = estimate) %>%
      arrange(adj_p)
  }

list_all_object_sizes <- 
  function() {
    
    
    
    objects(envir = .GlobalEnv) %>%
      enframe() %>%
      mutate(byte = sapply(value,
                           function(x) {
                             base::eval(parse(text = paste0("object.size(",
                                                            x,
                                                            ")")))
                           }),
             Mb = byte / 1e6,
             Gb = byte / 1e9) %>% 
      arrange(-Gb)
    
  }


color_mean <- 
  function(z) {
    col2rgb(z) %>%
      t() %>%
      as_tibble() %>%
      summarise(red = mean(red),
                green = mean(green),
                blue = mean(blue)) %$%
      rgb(red, green, blue, maxColorValue = 255)
    
    # rnorm(1)
    
    
  }

makeHexData <- function(x, y, z, bins, bin_fun) {
  
  h <- 
    hexbin::hexbin(x, y, bins, 
                   xbnds = range(x),
                   ybnds = range(y), 
                   IDs = TRUE)
  
  hexbin::hcell2xy(h) %>% 
    as_tibble() %>% 
    mutate(z = tapply(z, h@cID, FUN = bin_fun),
           cid = h@cell)
}


create_folder_structure <-
  function(dataset_id, dataset_run_id, main_folder, folders) {
    dir.create(main_folder, showWarnings = FALSE)
    dataset_folder <-
      tibble(dataset_id,
             dataset_run_id,
             dataset_path = paste0(main_folder, "/",
                                   dataset_id, "_",
                                   dataset_run_id))
    for(path in dataset_folder$dataset_path) dir.create(path, showWarnings = FALSE)
    folders_df <-
      enframe(folders,
              "name",
              "path") %>%
      expand_grid(dataset_id) %>%
      left_join(tibble(dataset_id,
                       dataset_run_id),
                by = "dataset_id") %>%
      left_join(dataset_folder,
                by = c("dataset_id",
                       "dataset_run_id")) %>%
      mutate(full_path = paste(dataset_path, path, sep = "/"))
    for(path in folders_df$full_path) dir.create(path, showWarnings = FALSE)
    out_paths <-
      set_names(dataset_folder$dataset_id,
                dataset_folder$dataset_id) %>%
      lapply(function(x) {
        clustering_info <-
          dataset_folder %>%
          filter(dataset_id == x)
        temp <-
          folders_df %>%
          filter(dataset_id == x)
        list(main = clustering_info$dataset_path,
             run_id = clustering_info$dataset_run_id) %>%
          append(split(temp$full_path,
                       temp$name))
      })
    return(out_paths)
  }

spread_colors <- 
  function(color, n, colorrange = 0.5) {
    
    if(n == 1) {
      return(color)
    } else {
      colorRamp(c("white", color, "black"))(seq(0.5 - colorrange / 2,
                                                0.5 + colorrange / 2, 
                                                length.out = n)) %>% 
        as_tibble() %$% 
        rgb(V1, V2, V3, maxColorValue = 255) 
      
    }
    
  }


complete_palette <- 
  function(pal_terms, pal, default_pal) {
    pal_terms %>% 
      enframe("i", "name") %>% 
      left_join(pal %>% 
                  enframe("name", "color"),
                by = c("name")) %>% 
      group_by(group = is.na(color)) %>% 
      mutate(color2 = default_pal(n_distinct(name))[row_number()]) %>% 
      ungroup() %>% 
      mutate(color = ifelse(is.na(color),
                            color2, 
                            color)) %>% 
      select(name, color) %>% 
      deframe()
    
  }


complete_palette_neighbor <- 
  function(pal_terms, pal, force_colors = "") {
    
    pal_terms %>% 
      enframe("i", "name") %>% 
      select(-i) %>% 
      expand_grid(term = names(pal)) %>% 
      distinct() %>% 
      left_join(pal %>% 
                  enframe("term", "color"),
                by = c("term")) %>% 
      distinct() %>% 
      mutate(strsim = stringdist::stringsim(name, term)) %>% 
      arrange(-strsim) %>% 
      group_by(name) %>% 
      slice(1) %>% 
      ungroup() %>% 
      select(name, color) %>% 
      mutate(color = case_when(name %in% names(force_colors) ~ force_colors[match(name, names(force_colors))],
                               T ~ color)) %>% 
      deframe()
    
  }

calculate_tau_score <- 
  function(wide_data) {
    max_exp <- 
      apply(wide_data,
            MARGIN = 1,
            function(x) max(x, na.rm = T))
    
    N <- 
      apply(wide_data,
            MARGIN = 1,
            function(x) length(which(!is.na(x))))
    
    expression_sum <- 
      wide_data %>% 
      sweep(MARGIN = 1, 
            STATS = max_exp, 
            FUN = `/`) %>% 
      {1 - .} %>% 
      apply(MARGIN = 1,
            function(x) sum(x, na.rm = T))
    
    
    tau_score <- 
      (expression_sum / (N - 1)) %>% 
      enframe("gene", "tau_score")
    
    tau_score
  }



ramp_color_bias <-
  function(color1, color2, biases) {
    require(tidyverse)
    colors_rgb <-
      colorRamp(c(color1, color2))(biases)
    rgb(colors_rgb[,1],
        colors_rgb[,2],
        colors_rgb[,3],
        maxColorValue = 255)
  }


cluster_long_data <-  
  function(long_data,
           distance_method = "euclidean",
           clustering_method = "ward.D2", 
           cluster_rows = T,
           cluster_cols = T) {
    suppressMessages(require(tidyverse))
    
    wide_data <- 
      long_data %>% 
      select(1:3) %>% 
      spread(2, 3) %>% 
      column_to_rownames(names(long_data)[1])
    
    order_row <- 
      rownames(wide_data)
    order_col <- 
      colnames(wide_data)
    
    
    if(cluster_rows) {
      order1 <- 
        wide_data %>% 
        dist(method = distance_method) %>% 
        hclust(method = clustering_method) %>% 
        with(labels[order])
    }
    
    if(cluster_cols) {
      order2 <- 
        wide_data %>% 
        t() %>% 
        dist(method = distance_method) %>% 
        hclust(method = clustering_method) %>% 
        with(labels[order])
    }
    
    long_data %>% 
      rename(v1 = 1, 
             v2 = 2,
             val = 3) %>% 
      mutate(v1 = factor(v1, order1),
             v2 = factor(v2, order2)) %>% 
      set_names(names(long_data))
    
  }
