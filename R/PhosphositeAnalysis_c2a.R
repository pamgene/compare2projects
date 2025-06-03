read_psites_from_folder <- function(folder, foldernum){
  files <- list.files(path = folder, pattern = "MTvC_|TT_", full.names = TRUE)
  lfcs <- files[grepl("_LogFC", files)]
  ps <- files[grepl("_p", files)]
  tercen <- files[grepl(".csv", files)]
  if (length(lfcs)){
    assign(paste0('lfcs', foldernum), lfcs, envir = parent.frame())
    assign(paste0('ps', foldernum), ps, envir = parent.frame())
  }
  if (length(tercen)){
    assign(paste0('tercen', foldernum), tercen, envir = parent.frame()) 
  }
}

are_lfcs_ps_paired <- function(lfcs, ps, foldernum){
  stopifnot("There are no files!" = length(lfcs) != 0)
  stopifnot("An LFC or p file is missing!" = length(lfcs) == length(ps))
  lfc_names <- sub("_LogFC.txt.*", "", lfcs)
  p_names <- sub("_p.txt.*", "", ps)
  stopifnot("LFC and p files are not paired!" = all(lfc_names == p_names))
}

is_data_in_2folders_paired_bn <- function(lfcs1, lfcs2, ps1, ps2){
  stopifnot("Number of files in the two folders are not equal!" = length(lfcs1) == length(lfcs2))
  lfc1_names <- sub(".*\\/(.+)_LogFC.txt.*", "\\1", lfcs1)
  lfc2_names <- sub(".*\\/(.+)_LogFC.txt.*", "\\1", lfcs2)
  p1_names <- sub(".*\\/(.+)_p.txt.*", "\\1", ps1)
  p2_names <- sub(".*\\/(.+)_p.txt.*", "\\1", ps2)
  stopifnot("LFC files in the two folders are not paired!" = all(lfc1_names == lfc2_names))
  stopifnot("p files in the two folders are not paired!" = all(p1_names == p2_names))
}

is_data_in_2folders_paired_tercen <- function(tercen1, tercen2){
  stopifnot("Number of files in the two folders are not equal!" = length(tercen1) == length(tercen2))
  tercen1_names <- sub(".*\\/(.+).csv.*", "\\1", tercen1)
  tercen2_names <- sub(".*\\/(.+).csv.*", "\\1", tercen2)
  stopifnot("Tercen files in the two folders are not paired!" = all(tercen1_names == tercen2_names))
}

clean_tercen_columns <- function(df) {
  colnames(df) <- sapply(str_split(colnames(df), "\\."), tail, 1)
  df
}

extract_phosphosite_data_mtvc_bionav <- function(lfc_file, p_file, pcutoff) {
  titles <- read_lines(p_file, n_max = 2)
  split_titles <- unlist(str_split(titles[2], "\t"))
  p_names <- c("c_cluster", "ID", "Uniprot", "Empty", split_titles[5:length(split_titles)])
  p_names <- p_names[p_names != ""]
  logfc <- read_delim(lfc_file, skip = 3, col_names = p_names, show_col_types = FALSE, delim = '\t')
  pvalues <- read_delim(p_file, skip = 3, col_names = p_names, show_col_types = FALSE, delim = '\t')
  ctrl_title <- split_titles[which(sapply(split_titles[5:length(split_titles)], function(title) {
    all(is.na(as.numeric(logfc[[title]])))
  }))[1] + 4]
  dfs <- lapply(split_titles[5:length(split_titles)], function(title) {
    if (title != "" && title != ctrl_title) {
      df <- logfc %>% select(ID, title)
      df$pvalue <- pvalues[[title]]
      colnames(df) <- c("ID", "LogFC", "P")
      df$LogFC <- as.numeric(df$LogFC)
      df$P <- as.numeric(df$P)
      df$Comparison <- paste(title, "vs", ctrl_title)
      df
    }
  })
  dfs <- bind_rows(Filter(Negate(is.null), dfs))
  dfs <- dfs[!str_detect(dfs$ID, "ART_0"),] %>% filter(P < pcutoff)
  dfs$Comparison <- as.factor(dfs$Comparison)
  dfs
}

extract_phosphosite_data_tt_bionav <- function(lfc_file, p_file, pcutoff, comp_from_fname) {
  titles <- read_lines(p_file, n_max = 2)
  split_titles <- unlist(str_split(titles[2], "\t"))
  if (length(split_titles) == 1) {
    lfc_cols <- c("ID", "Uniprot", "Sequence", "Empty", "LogFC")
    p_cols <- c("ID", "Uniprot", "Sequence", "Empty", "P")
    logfc <- read_delim(lfc_file, skip = 3, col_names = lfc_cols, show_col_types = FALSE, delim = '\t')
    p <- read_delim(p_file, skip = 3, col_names = p_cols, show_col_types = FALSE, delim = '\t')
    df <- logfc %>% select(ID, LogFC)
    df$P <- p$P
    df <- df %>% mutate(LogFC = as.numeric(LogFC), P = as.numeric(P), Comparison = comp_from_fname)
  } else {
    p_names <- c("ID", "Uniprot", "Sequence", "Empty", split_titles[5:length(split_titles)])
    p_names <- p_names[p_names != ""]
    logfc <- read_delim(lfc_file, skip = 3, col_names = p_names, show_col_types = FALSE, delim = "\t")
    pvalues <- read_delim(p_file, skip = 3, col_names = p_names, show_col_types = FALSE, delim = "\t")
    df <- lapply(split_titles[5:length(split_titles)], function(title) {
      if (title != "") {
        comp_df <- logfc %>% select(ID, title)
        comp_df$pvalue <- pvalues[[title]]
        colnames(comp_df) <- c("ID", "LogFC", "P")
        comp_df$LogFC <- as.numeric(comp_df$LogFC)
        comp_df$P <- as.numeric(comp_df$P)
        comp_df$Comparison <- title
        comp_df
      }
    })
    df <- bind_rows(Filter(Negate(is.null), df))
    df$Comparison <- as.factor(df$Comparison)
  }
  df <- df[!str_detect(df$ID, "ART_0"),] %>% filter(P < pcutoff)
  df
}

extract_phosphosite_data_mtvc_tercen <- function(filepath, pcutoff) {
  df <- read_delim(filepath, show_col_types = FALSE) %>% clean_tercen_columns()
  df <- df %>%
    rename(Comparison = 1) %>%
    mutate(variable = case_when(grepl("\\.LogFC", variable) ~ "LogFC", grepl("\\.pvalue", variable) ~ "P"))
  ctrl <- df %>% filter(is.na(value)) %>% distinct(Comparison) %>% pull()
  df %>% filter(Comparison != ctrl) %>%
    mutate(Comparison = paste(Comparison, "vs", ctrl), Comparison = as.factor(Comparison)) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    filter(P < pcutoff)
}

extract_phosphosite_data_tt_tercen <- function(filepath, pcutoff, comp_from_name) {
  df <- read_delim(filepath, show_col_types = FALSE) %>% clean_tercen_columns()
  if (length(colnames(df)) == 5) {
    df %>% rename(Comparison = 1, LogFC = delta, P = p) %>%
      mutate(Comparison = as.factor(Comparison)) %>%
      filter(P < pcutoff)
  } else {
    df %>% rename(LogFC = delta, P = p) %>%
      mutate(Comparison = comp_from_name) %>%
      filter(P < pcutoff)
  }
}

get_commonids <- function(sign1, sign2, fname, stat_type){
  conditions <- str_sort(as.character(unique(sign1$Comparison)))
  get_ids_per_condition <- function(sign){
    group <- aggregate(ID ~ Comparison, sign, paste0, collapse = " ")
    group1 <- setNames(lapply(group$ID, function(x) unlist(str_split(x, " "))), conditions)
    group1
  }
  sign1_ids <- get_ids_per_condition(sign1)
  sign2_ids <- get_ids_per_condition(sign2)
  common_ids <- map2(sign1_ids, sign2_ids, intersect)
  if (any(sapply(common_ids, length) == 0)){
    print(paste0("WARNING: ", names(which(sapply(common_ids, length) == 0)), ' conditions has no common ids.'))
  }
  n_common_ids_mat <- sapply(common_ids, length)
  n_sign1_ids_mat <- sapply(sign1_ids, length)
  n_sign2_ids_mat <- sapply(sign2_ids, length)
  df_123 <- data.frame(
    Comparison = names(n_common_ids_mat),
    n_common_peptides = n_common_ids_mat,
    data1_n_peptides = n_sign1_ids_mat,
    data2_n_peptides = n_sign2_ids_mat
  )
  write_csv(df_123, paste0("./results/result_common_peptides_", fname, "_", stat_type, "_numbers.csv"))
  common_ids
}

plot_common_data <- function(sign1, sign2, commonids, data1_name, data2_name, fname, stat_type){
  conditions <- str_sort(as.character(unique(sign1$Comparison)))
  conditions_data2 <- str_sort(as.character(unique(sign2$Comparison)))
  if (!identical(conditions, conditions_data2)){
    print("Conditions in data 1:"); print(conditions)
    print("Conditions in data 2:"); print(conditions_data2)
    stop(paste("Conditions in ", fname, " file do not match!\n",
               "conditions in data1: ", paste(conditions, collapse = ", "),
               "\nconditions in data2: ", paste(conditions_data2, collapse = ", "), sep = ""))
  }
  dflist <- lapply(seq_along(conditions), function(i){
    cond_df1 <- sign1[sign1$ID %in% commonids[[i]] & sign1$Comparison == conditions[i], c("ID", "Comparison", "LogFC")] %>%
      mutate(df = 'data1_LogFC') %>% pivot_wider(names_from = df, values_from = LogFC)
    cond_df2 <- sign2[sign2$ID %in% commonids[[i]] & sign2$Comparison == conditions[i], c("ID", "Comparison", "LogFC")] %>%
      mutate(df = 'data2_LogFC') %>% pivot_wider(names_from = df, values_from = LogFC)
    cond_df <- left_join(cond_df1, cond_df2, by = c("ID", "Comparison"))
    if (nrow(cond_df) != 0) cond_df else NULL
  })
  dflist <- Filter(Negate(is.null), dflist)
  if (length(dflist) > 0){
    commondf <- bind_rows(dflist)
    write_csv(commondf, paste0("./results/result_common_peptides_", fname, "_", stat_type, ".csv"))
    title <- paste("Significant peptides in common - ", fname)
    filename <- paste0("./results/result_common_peptides_", fname, "_", stat_type, ".png")
    w <- ifelse(length(conditions) == 1, 10, 16)
    h <- 6
    wrapper <- function(label, dev_width = dev.size("in")[1], dev_scaler = 5) paste(strwrap(label, dev_width * dev_scaler), collapse = "\n")
    print(paste("Plotting ", fname,"...", sep = ""))
    myplot <- ggplot(commondf, aes(data1_LogFC, data2_LogFC)) +
      geom_point() +
      facet_wrap(~ Comparison) +
      ylab(data2_name) +
      xlab(data1_name) +
      xlim(-2, 2) +
      ylim(-2,2) +
      ggtitle(wrapper(title))
    ggsave(filename, myplot, width = w, height = h, units = "cm")
  }
}

get_peptide_results_bn <- function(lfcs1, ps1, lfcs2, ps2, pcutoff, data1_name = "data1", data2_name = "data2"){
  for (i in seq_along(lfcs1)){
    comparison <- sub(".*\\/(.+)_LogFC.txt.*", "\\1", lfcs1[[i]])
    print(paste("Working on file:", comparison))
    stat_type <- ifelse(grepl("MTvC", lfcs1[[i]]), "MTvC", "TT")
    sign1 <- if (stat_type == "MTvC") extract_phosphosite_data_mtvc_bionav(lfcs1[[i]], ps1[[i]], pcutoff)
             else extract_phosphosite_data_tt_bionav(lfcs1[[i]], ps1[[i]], pcutoff, comparison)
    sign2 <- if (stat_type == "MTvC") extract_phosphosite_data_mtvc_bionav(lfcs2[[i]], ps2[[i]], pcutoff)
             else extract_phosphosite_data_tt_bionav(lfcs2[[i]], ps2[[i]], pcutoff, comparison)
    if (nrow(sign1) == 0 | nrow(sign2) == 0){
      print("WARNING: There are no significant phosphosites in at least one data")
    } else {
      commonids <- try(get_commonids(sign1, sign2, comparison, stat_type))
      try(plot_common_data(sign1, sign2, commonids, data1_name, data2_name, comparison, stat_type))
    }
  }
}

get_peptide_results_tercen <- function(tercen1, tercen2, pcutoff, data1_name = "data1", data2_name = "data2"){
  for (i in seq_along(tercen1)){
    comparison <- sub(".*\\/(.+).csv.*", "\\1", tercen1[[i]])
    print(paste("Working on file:", comparison))
    stat_type <- ifelse(grepl("MTvC", tercen1[[i]]), "MTvC", "TT")
    sign1 <- if (stat_type == "MTvC") extract_phosphosite_data_mtvc_tercen(tercen1[[i]], pcutoff)
             else extract_phosphosite_data_tt_tercen(tercen1[[i]], pcutoff, comparison)
    sign2 <- if (stat_type == "MTvC") extract_phosphosite_data_mtvc_tercen(tercen2[[i]], pcutoff)
             else extract_phosphosite_data_tt_tercen(tercen2[[i]], pcutoff, comparison)
    if (nrow(sign1) == 0 | nrow(sign2) == 0){
      print("WARNING: There are no significant phosphosites in at least one data")
    } else {
      commonids <- try(get_commonids(sign1, sign2, comparison, stat_type))
      try(plot_common_data(sign1, sign2, commonids, data1_name, data2_name, comparison, stat_type))
    }
  }
}



