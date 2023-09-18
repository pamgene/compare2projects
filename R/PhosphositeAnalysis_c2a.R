
read_psites_from_folder <- function(folder, foldernum){
  files <- list.files(path = folder, pattern = "MTvC_|TT_", full.names = TRUE)
  lfcs <- files[grepl("_LogFC", files)]
  ps <- files[grepl("_p", files)]
  tercen <- files[grepl(".csv", files)]
  if (!is_empty(lfcs)){
    assign(paste0('lfcs', foldernum), lfcs, envir = parent.frame())
    assign(paste0('ps', foldernum), ps, envir = parent.frame())
  }
  if (!is_empty(tercen)){
    assign(paste0('tercen', foldernum), tercen, envir = parent.frame()) 
  }
}

are_lfcs_ps_paired <- function(lfcs, ps, foldernum){
  # 1. stop if an lfc or p is missing 
  stopifnot("There are no files!" = length(lfcs) != 0)
  stopifnot("An LFC or p file is missing!" = length(lfcs) == length(ps))
  
  # 2. stop if the names don't match
  lfc_names = c()
  p_names = c()
  for (i in seq_along(lfcs)){
    #get name
    lfc_name <- sub("_LogFC.txt.*", "", lfcs[[i]])
    p_name <- sub("_p.txt.*", "", ps[[i]])
    lfc_names <- append(lfc_names, lfc_name)
    p_names <- append(p_names, p_name)
    
  }
  stopifnot("LFC and p files are not paired!" = all(lfc_names == p_names))
  #stopifnot(paste0("LFC and p files in folder ", foldernum, " are not paired!") = all(lfc_names == p_names))
}

is_data_in_2folders_paired_bn <- function(lfcs1, lfcs2, ps1, ps2){
  # 1. stop if number of lfcs are not the same
  stopifnot("Number of files in the two folders are not equal!" = length(lfcs1) == length(lfcs2))
  
  # 2. check if names are the same
  
  lfc1_names = c()
  lfc2_names = c()
  p1_names = c()
  p2_names = c()
  for (i in seq_along(lfcs1)){
    #get name
    lfc1_name <- sub(".*\\/(.+)_LogFC.txt.*", "\\1", lfcs1[[i]])
    lfc2_name <- sub(".*\\/(.+)_LogFC.txt.*", "\\1", lfcs2[[i]])
    p1_name <- sub(".*\\/(.+)_p.txt.*", "\\1", ps1[[i]])
    p2_name <- sub(".*\\/(.+)_p.txt.*", "\\1", ps2[[i]])
    lfc1_names <- append(lfc1_names, lfc1_name)
    lfc2_names <- append(lfc2_names, lfc2_name)
    p1_names <- append(p1_names, p1_name)
    p2_names <- append(p2_names, p2_name)
  }
  
  stopifnot("LFC files in the two folders are not paired!" = all(lfc1_names == lfc2_names))
  stopifnot("p files in the two folders are not paired!" = all(p1_names == p2_names))
}

is_data_in_2folders_paired_tercen <- function(tercen1, tercen2){
  # 1. stop if number of files are not the same
  stopifnot("Number of files in the two folders are not equal!" = length(tercen1) == length(tercen2))
  
  # 2. check if names are the same
  tercen1_names = c()
  tercen2_names = c()
  
  for (i in seq_along(tercen1)){
    #get name
    tercen1_name <- sub(".*\\/(.+).csv.*", "\\1", tercen1[[i]])
    tercen2_name <- sub(".*\\/(.+).csv.*", "\\1", tercen2[[i]])
    tercen1_names <- append(tercen1_names, tercen1_name)
    tercen2_names <- append(tercen2_names, tercen2_name)
  }
  stopifnot("Tercen files in the two folders are not paired!" = all(tercen1_names == tercen2_names))
}


clean_tercen_columns <- function(df) {
  cols <- colnames(df)
  split <- str_split(cols, pattern = "\\.")
  cols <- sapply(split, tail, 1)
  colnames(df) <- cols
  return(df)
}



extract_phosphosite_data_mtvc_bionav <- function(lfc_file, p_file, pcutoff) {
  #print(lfc_file)
  #print(p_file)
  titles <- read_lines(p_file, n_max = 2)
  #print(paste("titles", titles))
  split_titles <- unlist(str_split(titles[2], "\t"))
  #print(paste("split_titles", split_titles))
  p_names <- c("c_cluster", "ID", "Uniprot", "Empty", split_titles[5:length(split_titles)])
  p_names <- p_names[p_names != ""]
  
  logfc <- read_delim(lfc_file, skip = 3, col_names = p_names, show_col_types = FALSE, delim = '\t')
  pvalues <- read_delim(p_file, skip = 3, col_names = p_names, show_col_types = FALSE, delim = '\t')
  
  dfs <- list()
  ctrl_title <- ""
  
  # first find which column is control, e.g., all NAs
  for (title in split_titles[5:length(split_titles)]) {
    if (title != "") {
      df <- logfc %>% select(ID, title)
      colnames(df) <- c("ID", "LogFC")
      df$LogFC <- as.numeric(df$LogFC)
      
      if (suppressWarnings(is.na(any(df$LogFC)))) {
        ctrl_title <- title
      }
    }
  }
  
  for (title in split_titles[5:length(split_titles)]) {
    if ((title != ctrl_title) & (title != "")) {
      df <- logfc %>% select(ID, title)
      df$pvalue <- pvalues %>%
        select(title) %>%
        pull()
      colnames(df) <- c("ID", "LogFC", "P")
      df$LogFC <- as.numeric(df$LogFC)
      df$P <- as.numeric(df$P)
      
      comparison <- paste(title, "vs", ctrl_title)
      df$Comparison <- comparison
      dfs[[length(dfs) + 1]] <- df
    }
  }
  
  dfs <- bind_rows(dfs)
  dfs$Comparison <- as.factor(dfs$Comparison)
  
  dfs <- dfs[!str_detect(dfs$ID, "ART_0"),]
  
  dfs_s <- dfs %>% filter(P<pcutoff)
  
  return(dfs_s)
}


extract_phosphosite_data_tt_bionav <- function(lfc_file, p_file, pcutoff, comp_from_fname) {
  
  # identify whether it's a supergroup file 
  titles <- read_lines(p_file, n_max = 2)
  split_titles <- unlist(str_split(titles[2], "\t"))
  
  if (length(split_titles) == 1) {
    # not supergroup 
    lfc_cols <- c("ID", "Uniprot", "Sequence", "Empty", "LogFC")
    p_cols <- c("ID", "Uniprot", "Sequence", "Empty", "P")
    
    logfc <- read_delim(lfc_file, skip = 3, col_names = lfc_cols, show_col_types = FALSE, delim = '\t')
    p <- read_delim(p_file, skip = 3, col_names = p_cols, show_col_types = FALSE, delim = '\t')
    df <- logfc %>% select(ID, LogFC)
    df$P <- p %>%
      select(P) %>%
      pull()
    df <- df %>% mutate(
      LogFC = as.numeric(LogFC),
      P = as.numeric(P)
    )
    df$Comparison <- comp_from_fname
    
  } else {
    # Supergroup
    p_names <- c("ID", "Uniprot", "Sequence", "Empty", split_titles[5:length(split_titles)])
    p_names <- p_names[p_names != ""]

    logfc <- read_delim(lfc_file, skip = 3, col_names = p_names, show_col_types = FALSE, delim = "\t")
    pvalues <- read_delim(p_file, skip = 3, col_names = p_names, show_col_types = FALSE, delim = "\t")
    
    df <- list()
    for (title in split_titles[5:length(split_titles)]) {
      if ((title != "")) {
        comp_df <- logfc %>% select(ID, title)
        
        comp_df$pvalue <- pvalues %>%
          select(title) %>%
          pull()
        
        colnames(comp_df) <- c("ID", "LogFC", "P")
        comp_df$LogFC <- as.numeric(comp_df$LogFC)
        comp_df$P <- as.numeric(comp_df$P)
        
        comparison <- paste(title)
        comp_df$Comparison <- comparison
        df[[length(df) + 1]] <- comp_df
      }
    }
    df <- bind_rows(df)
    df$Comparison <- as.factor(df$Comparison)
  }
  
  
  df <- df[!str_detect(df$ID, "ART_0"),]
  df_s <- df %>% filter(P<pcutoff)
  
  return(df_s)
  
}

extract_phosphosite_data_mtvc_tercen <- function(filepath, pcutoff) {
  
  df <- read_delim(filepath, show_col_types = FALSE)
  df <- clean_tercen_columns(df)
  df <- df %>%
    rename(Comparison = 1) %>%
    mutate(variable = case_when(grepl("\\.LogFC", variable) ~ "LogFC", grepl("\\.pvalue", variable) ~ "P"))
  ctrl <- df %>%
    filter(is.na(value)) %>%
    distinct(Comparison) %>%
    pull()
  
  df <- df %>% filter(Comparison != ctrl) %>%
    mutate(Comparison = paste(Comparison, "vs", ctrl),
           Comparison = as.factor(Comparison)) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    filter(P<pcutoff)
  return(df)
}


extract_phosphosite_data_tt_tercen <- function(filepath, pcutoff, comp_from_name) {
  df <- read_delim(filepath, show_col_types = FALSE)
  df <- clean_tercen_columns(df)
  
  if (length(colnames(df)) == 5) {
    # supergroup
    df <- df %>%
      rename(Comparison = 1, LogFC = delta, P = p) %>%
      mutate(Comparison = as.factor(Comparison)) %>%
      filter(P<pcutoff)
  } else {
    # no supergroup
    df <- df %>%
      rename(LogFC = delta, P = p) %>%
      mutate(Comparison = comp_from_name) %>%
      filter(P<pcutoff) 
  }

  return(df)
}



get_commonids <- function(sign1, sign2, fname, stat_type){
  
  conditions <- as.character(unique(sign1$Comparison))
  conditions <- str_sort(conditions)
  
  get_ids_per_condition <- function(sign){
    group <- aggregate(ID~Comparison, sign, paste0, collapse = " ")
    group1 <- sapply(group$ID, str_split, " ")
    names(group1) <- conditions
    return(group1)
  }
  
  sign1_ids <- get_ids_per_condition(sign1)
  sign2_ids <- get_ids_per_condition(sign2)
  common_ids <- map2(sign1_ids, sign2_ids, intersect)
  
  if (any(lapply(common_ids, length) == 0)){
    print(paste0("WARNING: ", names(lapply(common_ids, length) == 0), ' conditions has no common ids.'))
  }
  
  n_common_ids_mat <- sapply(common_ids, length)
  n_sign1_ids_mat <- sapply(sign1_ids, length)
  n_sign2_ids_mat <- sapply(sign2_ids, length)
  df_common <- data.frame(n_common_peptides = n_common_ids_mat) %>%
    rownames_to_column("Comparison")
  df_sign1 <- data.frame(data1_n_peptides = n_sign1_ids_mat) %>%
    rownames_to_column("Comparison")
  df_sign2 <- data.frame(data2_n_peptides = n_sign2_ids_mat) %>%
    rownames_to_column("Comparison")
  
  df_12 <- left_join(df_common, df_sign1, by = "Comparison")
  df_123 <- left_join(df_12, df_sign2, by = "Comparison")
  
  write_csv(df_123, paste("./results/result_common_peptides_", fname, "_", stat_type, "_numbers", ".csv", sep = ""))
  
  return(common_ids)
  
}


plot_common_data <- function(sign1, sign2, commonids, data1_name, data2_name, fname, stat_type){
  conditions <- as.character(unique(sign1$Comparison)) %>% str_sort()
  conditions_data2 <- as.character(unique(sign2$Comparison)) %>% str_sort()
  if (any(conditions != conditions_data2)){
    print("Conditions in data 1:")
    print(conditions)
    print("Conditions in data 2:")
    print(conditions_data2)
    stop(paste("Conditions in ", fname, " file do not match!\n", 
               "conditions in data1: ", paste(conditions, collapse = ", "),
               "\nconditions in data2: ", paste(conditions_data2, collapse = ", "), sep = ""))
  } 
  
  
  dflist <- as.list(rep(NA, length(conditions)))

  # for each condition (conditions is a vector)
  # commonids is a list with as many elements as conditions
  for (i in seq_along(conditions)){
    cond_df1 <- sign1[sign1$ID %in% commonids[[i]], c("ID", "Comparison", "LogFC")] %>%
      filter(Comparison == conditions[i]) %>%
      mutate(df = 'data1_LogFC')%>% 
      pivot_wider(names_from = df, values_from = LogFC)
    cond_df2 <- sign2[sign2$ID %in% commonids[[i]], c("ID", "Comparison", "LogFC")] %>%
      filter(Comparison == conditions[i]) %>%
      mutate(df = 'data2_LogFC')%>% 
      pivot_wider(names_from = df, values_from = LogFC)
    cond_df <- left_join(cond_df1, cond_df2, by = c("ID", "Comparison"))
    if (nrow(cond_df)!=0){
      dflist[[i]] <- cond_df
    }
    
  }

  dflist <- dflist[!is.na(dflist)]
  if (length(dflist)>0){
    commondf <- dflist %>% purrr::reduce(rbind)
    write_csv(commondf, paste("./results/result_common_peptides_", fname, "_", stat_type, ".csv", sep = ""))
    
    # params for plotting
    title <- paste("Significant peptides in common - ", fname)
    filename <- paste("./results/result_common_peptides_", fname, "_", stat_type, ".png", sep = '')
    
    w <- ifelse(length(conditions) == 1, 10, 16) 
    # since TT is usually 1. 
    # but if there is only 1 condition plotted from the two, it should also be one (think of solution!)
    #w <- 12 
    h <- 6
    
    wrapper <- function(label, dev_width = dev.size("in")[1], dev_scaler = 5)  {   
      paste(strwrap(label, dev_width * dev_scaler), collapse = "\n") 
    }
    
    print(paste("Plotting ", fname,"...", sep = ""))
    
    myplot <- ggplot(commondf, aes(data1_LogFC, data2_LogFC)) +
      geom_point() +
      facet_wrap(~ Comparison) +
      ylab(data2_name) +
      xlab(data1_name) +
      xlim(-2, 2) +
      ylim(-2,2) +
      ggtitle(wrapper(title))
    
    ggsave(filename, myplot,
           width = w, height = h, units = "cm")
    
  }
  
}



get_peptide_results_bn <- function(lfcs1, ps1, lfcs2, ps2, pcutoff,
                                data1_name = "data1", data2_name = "data2"){
  
  for (i in seq_along(lfcs1)){
    # get comparison
    comparison <- sub(".*\\/(.+)_LogFC.txt.*", "\\1", lfcs1[[i]])
    print(paste("Working on file:", comparison))
    
    # get assay type
    assay_type <- ifelse(grepl("_PTK_", lfcs1[[i]]), "PTK", "STK")
    
    # get stat type
    stat_type <- ifelse(grepl("MTvC", lfcs1[[i]]), "MTvC", "TT")
    
    
    if (stat_type == "MTvC"){
      sign1 <- extract_phosphosite_data_mtvc_bionav(lfcs1[[i]], ps1[[i]], pcutoff)
      sign2 <- extract_phosphosite_data_mtvc_bionav(lfcs2[[i]], ps2[[i]], pcutoff)
    } else if (stat_type == "TT"){
      sign1 <- extract_phosphosite_data_tt_bionav(lfcs1[[i]], ps1[[i]], pcutoff, comparison)
      sign2 <- extract_phosphosite_data_tt_bionav(lfcs2[[i]], ps2[[i]], pcutoff, comparison)
    }

    if (dim(sign1)[1] == 0 | dim(sign2)[1] == 0){
      print("WARNING: There are no significant phosphosites in at least one data")
    }  else {
      commonids <- try(get_commonids(sign1, sign2, comparison, stat_type))
      try(plot_common_data(sign1, sign2, commonids, data1_name, data2_name, comparison, stat_type))
    }
    
    # commonids <- get_commonids(sign1, sign2, fname, stat_type)
    # print("commonids")
    # print(commonids)
    # plot_common_data(sign1, sign2, commonids, data1_name, data2_name, fname, stat_type)
    
  }
}

get_peptide_results_tercen <- function(tercen1, tercen2, pcutoff,
                                   data1_name = "data1", data2_name = "data2"){
  
  for (i in seq_along(tercen1)){
    # get comparison
    comparison <- sub(".*\\/(.+).csv.*", "\\1", tercen1[[i]])
    print(paste("Working on file:", comparison))
    
    # get assay type
    assay_type <- ifelse(grepl("_PTK_", tercen1[[i]]), "PTK", "STK")
    
    # get stat type
    stat_type <- ifelse(grepl("MTvC", tercen1[[i]]), "MTvC", "TT")
    
    
    if (stat_type == "MTvC"){
      sign1 <- extract_phosphosite_data_mtvc_tercen(tercen1[[i]], pcutoff)
      sign2 <- extract_phosphosite_data_mtvc_tercen(tercen2[[i]], pcutoff)
    } else if (stat_type == "TT"){
      sign1 <- extract_phosphosite_data_tt_tercen(tercen1[[i]], pcutoff, comparison)
      sign2 <- extract_phosphosite_data_tt_tercen(tercen2[[i]], pcutoff, comparison)
    }
    if (dim(sign1)[1] == 0 | dim(sign2)[1] == 0){
      print("WARNING: There are no significant phosphosites in at least one data")
    }  else {
      commonids <- try(get_commonids(sign1, sign2, comparison, stat_type))
      try(plot_common_data(sign1, sign2, commonids, data1_name, data2_name, comparison, stat_type))
    }
  }
}



