

read_uka_from_folder <- function(folder){
  files <- list.files(path = path, pattern = "UKA_", 
                      full.names = TRUE)
  res <- lapply(files, read_delim, show_col_types = F)
  return(res)
  
}


parse_uka <- function(uka, fscorecutoff){
  long <- uka %>% 
    filter(`Median Final score` >= fscorecutoff) %>%
    select(`Kinase Name`, `Mean Specificity Score`, `Mean Significance Score`,
           `Median Final score`, `Mean Kinase Statistic`) 
  long <- long %>%
    pivot_longer(cols = colnames(long[,2:dim(long)[2]]), 
                 names_to = "scoretype", values_to = "value")
  
  return(long)
}


filter_for_common <- function(uka, commonids){
  commonuks <- uka %>% filter(`Kinase Name` %in% commonids)
  return(commonuks)
}




compare_two_uka <- function(ukapath1, ukapath2, data1_name, data2_name, fscorecutoff = 1.3, fname){
  uka1 <- read_delim(ukapath1, show_col_types = F)
  uka2 <- read_delim(ukapath2, show_col_types = F)
  parsed1 <- parse_uka(uka1, fscorecutoff)
  parsed2 <- parse_uka(uka2, fscorecutoff)
  
  ids1 <- unique(parsed1$"Kinase Name")
  ids2 <- unique(parsed2$"Kinase Name")
  commonids <- intersect(ids1, ids2)
  
  if (length(commonids) == 0){
    print("WARNING: There are no common kinases in the data")
  }
  
  uka1_common <- filter_for_common(parsed1, commonids) %>% dplyr::rename("data1" = "value")
  uka2_common <- filter_for_common(parsed2, commonids) %>% dplyr::rename("data2" = "value")
  
  commondf <- left_join(uka1_common, uka2_common)
  
  if (dim(commondf)[1] > 0){
    write_csv(commondf, paste("./result_common_kinases_", fname, ".csv", sep = ""))
  }
  
  wrapper <- function(label, dev_width = dev.size("in")[1], dev_scaler = 5)  {   
    paste(strwrap(label, dev_width * dev_scaler), collapse = "\n") 
  }
  
  if (dim(commondf)[1] > 0){
    print(paste("Plotting UKA: ", fname,"...", sep = ""))
    ukaplot <- ggplot(commondf) +
      geom_point(aes(x = data1, y = data2)) +
      geom_point(aes(x = data2, y = data1), alpha = 0) + # makes x and y axes the same scale
      scale_alpha(range = c(0, 1)) +
      facet_wrap(~ scoretype, scales = 'free') +
      ylab(data2_name) +
      xlab(data1_name) +
      ggtitle(wrapper(paste("Significant kinases in common - ", fname)))
    
    
    ggsave(paste("./result_common_kinases_", fname, ".png", sep = ''), ukaplot,
           width = 12, height = 10, units = "cm")
    
  }
  
  
  
}


get_uka_results <- function(uka_files1, uka_files2, data1_name, data2_name, fscorecutoff){
  for (i in seq_along(uka_files1)){
    # get comparison
    comparison <- sub(".*\\/(.+).(txt|csv).*", "\\1", uka_files1[[i]])
    print(paste("Working on UKA file:", comparison))
    # get assay type
    assay_type <- ifelse(grepl("_PTK_", uka_files1[[i]]), "PTK", "STK")
    
    compare_two_uka(ukapath1 = uka_files1[[i]], ukapath2 = uka_files2[[i]],
                    data1_name = data1_name, data2_name = data2_name,
                    fname = comparison)
  }
  
}




