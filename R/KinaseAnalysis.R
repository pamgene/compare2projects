read_uka_from_folder <- function(folder){
  files <- list.files(path = folder, pattern = "UKA_", full.names = TRUE)
  lapply(files, read_delim, show_col_types = FALSE)
}

filter_for_common <- function(uka, commonids){
  uka %>% filter(`Kinase Name` %in% commonids)
}

clean_tercen_columns <- function(df) {
  colnames(df) <- sapply(str_split(colnames(df), "\\."), tail, 1)
  df
}

# Helper for wrapping plot titles
wrapper <- function(label, dev_width = dev.size("in")[1], dev_scaler = 5)  {   
  paste(strwrap(label, dev_width * dev_scaler), collapse = "\n") 
}

# Helper for saving plots
save_uka_plot <- function(plot, filename, width = 12, height = 10, units = "cm") {
  ggsave(filename, plot, width = width, height = height, units = units)
}

# Helper to get common and unique kinases and summary row
get_common_and_unique_kinases <- function(parsed1, parsed2, comparison = NA, data1_name = NULL, data2_name = NULL) {
  ids1 <- unique(parsed1$"Kinase Name")
  ids2 <- unique(parsed2$"Kinase Name")
  commonids <- intersect(ids1, ids2)
  only1 <- setdiff(ids1, ids2)
  only2 <- setdiff(ids2, ids1)
  overlap_coef <- length(commonids) / min(length(ids1), length(ids2))
  size_ratio <- min(length(ids1), length(ids2)) / max(length(ids1), length(ids2))
  # Use setNames to dynamically assign column names
  if (!is.null(data1_name) && !is.null(data2_name)) {
    summary_row <- data.frame(
      comparison = comparison,
      overlap_coef = overlap_coef,
      size_ratio = size_ratio
    )
    summary_row[[paste0("n_only_", data1_name)]] <- length(only1)
    summary_row[[paste0("n_only_", data2_name)]] <- length(only2)
  } else {
    summary_row <- data.frame(
      comparison = comparison,
      overlap_coef = overlap_coef,
      size_ratio = size_ratio,
      n_only_data1 = length(only1),
      n_only_data2 = length(only2)
    )
  }
  list(commonids = commonids, only1 = only1, only2 = only2, summary_row = summary_row)
}

# Helper to collect commondf for plotting/writing
collect_commondf <- function(parsed1, parsed2, commonids, extra_cols = list()) {
  uka1_common <- filter_for_common(parsed1, commonids) %>% dplyr::rename("data1" = "value")
  uka2_common <- filter_for_common(parsed2, commonids) %>% dplyr::rename("data2" = "value")
  commondf <- left_join(uka1_common, uka2_common)
  for (colname in names(extra_cols)) {
    commondf[[colname]] <- extra_cols[[colname]]
  }
  commondf
}

parse_uka <- function(uka, fscorecutoff){
  uka %>%
    filter(`Median Final score` >= fscorecutoff) %>%
    select(`Kinase Name`, `Mean Specificity Score`, `Mean Significance Score`,
           `Median Final score`, `Median Kinase Statistic`) %>%
    pivot_longer(-`Kinase Name`, names_to = "scoretype", values_to = "value")
}

# Helper to parse and filter a UKA file for a given contrast (optional)
parse_and_filter_uka <- function(uka, fscorecutoff, contrast = NULL) {
  if (!is.null(contrast) && "Sgroup_contrast" %in% colnames(uka)) {
    # Try to coerce both to character for comparison
    uka <- uka %>% filter(as.character(Sgroup_contrast) == as.character(contrast))
  }
  result <- parse_uka(uka, fscorecutoff)
  # Print head to shiny interface
  message("Head of parsed UKA:")
  message(capture.output(print(head(result))))
  result
}

make_kin_stat_barplot <- function(result_with_scores, min_comps, data1_name, data2_name){
  robust_data_to_plot <- result_with_scores %>%
    filter(n_comparisons > min_comps, robust_change, scoretype != "Median Final score") %>%
    select(-robust_change) %>%
    mutate(Kinase_ncomp = paste0(`Kinase Name`, " - ", n_comparisons))

  # Label aliases for scoretype
  scoretype_labels <- c(
    "Mean Specificity Score" = "Spec Score",
    "Median Kinase Statistic" = "Kin Stat",
    "Mean Significance Score" = "Sig Score",
    "Median Final score" = "Final Score"
  )
  # Only use labels for present levels
  present_levels <- intersect(names(scoretype_labels), unique(robust_data_to_plot$scoretype))
  robust_data_to_plot$scoretype <- factor(
    robust_data_to_plot$scoretype,
    levels = present_levels,
    labels = scoretype_labels[present_levels]
  )

  # Dynamically determine plot size
  n_kinases <- length(unique(robust_data_to_plot$`Kinase Name`))
  n_scoretypes <- length(unique(robust_data_to_plot$scoretype))
  width <- max(16, n_scoretypes * 4)
  height <- max(10, n_kinases * 1.5)

  # Add robustness thresholds as caption
  min_comps <- 4
  mean_thresh <- 0.2
  sd_thresh <- 0.3
  threshold_caption <- paste0(
    "Kinase name - number of comparisons \nRobustness thresholds: min_comps = ", min_comps,
    ", |mean diff| > ", mean_thresh,
    ", SD < ", sd_thresh,
    '\n error bars = SD of mean difference'
  )

  # Remove ordering by spec score; use default order
  robust_data_to_plot$Kinase_ncomp <- factor(robust_data_to_plot$Kinase_ncomp)

  robust_kin_plot <- ggplot(robust_data_to_plot, aes(x = scoretype, y = mean_difference, fill = scoretype)) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_errorbar(
      aes(ymin = mean_difference - sd_difference, ymax = mean_difference + sd_difference),
      width = 0.2,
      position = position_dodge(width = 0.9)
    ) +
    facet_wrap(~ Kinase_ncomp, scales = "fixed", ncol = 5) +
    labs(
      title = paste0("Robustly changing kinases (", data2_name, " - ", data1_name, ")"),
      y = "Mean of Score Differences",
      x = "Score Type",
      caption = threshold_caption
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.caption = element_text(hjust = 0)
    )

  # Save the plot as PNG in results folder
  ggsave(
    filename = "./results/common_kinases_with_robust_change.png",
    plot = robust_kin_plot,
    width = width,
    height = height,
    units = "cm"
  )
}

make_per_kin_stats <- function(df, data1_name, data2_name){
  # 1. summarize mean and sd of difference by kinase and scoretype
  kinase_summary <- df %>%
    group_by(`Kinase Name`, scoretype) %>%
    summarize(
      n_comparisons = n(),
      mean_difference = mean(difference, na.rm = TRUE),
      sd_difference = sd(difference, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # drop_na(sd_difference) %>% # if sd is NA, there is only 1 comparison
    arrange(`Kinase Name`, scoretype, desc(mean_difference))
  write_csv(kinase_summary, "./results/common_kinases_difference_summary.csv")

  # Flag kinases as "robustly changing" based on mean and sd thresholds
  min_comps <- 4
  mean_thresh <- 0.2
  sd_thresh <- 0.3

  df_flagged <- kinase_summary %>%
    filter(n_comparisons >= min_comps) %>%
    mutate(robust_change = abs(mean_difference) > mean_thresh & sd_difference < sd_thresh)
  df_flagged_sum <- df_flagged %>%
    group_by(`Kinase Name`) %>%
    summarise(
      num_robust_scores = sum(robust_change),
      has_any_robust = any(robust_change),
      total_scores = n()) %>%
    filter(has_any_robust)
  result_with_scores <- df_flagged %>%
    semi_join(df_flagged_sum, by = "Kinase Name") %>%
    arrange(`Kinase Name`, sd_difference)

  write_csv(result_with_scores, "./results/common_kinases_flagged.csv")
  # 2. make barplot of mean difference by kinase and scoretype
  make_kin_stat_barplot(result_with_scores = result_with_scores, min_comps = min_comps, data1_name = data1_name, data2_name = data2_name)
}


compare_two_uka <- function(ukapath1, ukapath2, data1_name, data2_name, fscorecutoff, fname){
  uka1 <- read_delim(ukapath1, show_col_types = F) %>% clean_tercen_columns()
  uka2 <- read_delim(ukapath2, show_col_types = F) %>% clean_tercen_columns()
  parsed1 <- parse_and_filter_uka(uka1, fscorecutoff)
  parsed2 <- parse_and_filter_uka(uka2, fscorecutoff)
  
  kinases_info <- get_common_and_unique_kinases(parsed1, parsed2)
  commonids <- kinases_info$commonids
  
  if (length(commonids) == 0){
    print("WARNING: There are no common kinases in the data")
  }
  
  commondf <- collect_commondf(parsed1, parsed2, commonids)
  
  if (nrow(commondf) > 0){
    write_csv(commondf, paste("./results/common_kinases_", fname, ".csv", sep = ""))
    print(paste("Plotting UKA: ", fname,"...", sep = ""))
    ukaplot <- ggplot(commondf) +
      geom_point(aes(x = data1, y = data2)) +
      geom_point(aes(x = data2, y = data1), alpha = 0) + # makes x and y axes the same scale
      scale_alpha(range = c(0, 1)) +
      facet_wrap(~ scoretype, scales = 'free') +
      ylab(data2_name) +
      xlab(data1_name) +
      ggtitle(wrapper(paste("Significant kinases in common - ", fname)))
    save_uka_plot(ukaplot, paste("./results/common_kinases_", fname, ".png", sep = ''))
  }
}

process_all_vs_all_uka <- function(uka_file1, uka_file2, fscorecutoff, data1_name, data2_name) {
  uka1 <- read_delim(uka_file1, show_col_types = FALSE) %>% clean_tercen_columns()
  uka2 <- read_delim(uka_file2, show_col_types = FALSE) %>% clean_tercen_columns()
  if (!("Sgroup_contrast" %in% colnames(uka1)) || !("Sgroup_contrast" %in% colnames(uka2))) {
    warning("Sgroup_contrast column not found in one or both UKA files for all_vs_all mode.")
    return()
  }
  contrasts <- intersect(unique(uka1$Sgroup_contrast), unique(uka2$Sgroup_contrast))
  # Output contrasts to shiny interface
  message(sprintf("Contrasts for comparison %s: %s", comparison, paste(contrasts, collapse = ", ")))
  for (contrast in contrasts) {
    print(paste("Working on UKA contrast:", contrast))
    parsed1 <- parse_and_filter_uka(uka1, fscorecutoff, contrast)
    parsed2 <- parse_and_filter_uka(uka2, fscorecutoff, contrast)
    kinases_info <- get_common_and_unique_kinases(parsed1, parsed2, contrast = contrast)
    commonids <- kinases_info$commonids
    if (length(commonids) == 0){
      print("WARNING: There are no common kinases in the data")
    }
    commondf <- collect_commondf(parsed1, parsed2, commonids)
    if (nrow(commondf) > 0){
      outname <- paste0("./results/result_common_kinases_", contrast, ".csv")
      write_csv(commondf, outname)
      print(paste("Plotting UKA: ", contrast,"...", sep = ""))
      ukaplot <- ggplot(commondf) +
        geom_point(aes(x = data1, y = data2)) +
        geom_point(aes(x = data2, y = data1), alpha = 0) +
        scale_alpha(range = c(0, 1)) +
        facet_wrap(~ scoretype, scales = 'free') +
        ylab(data2_name) +
        xlab(data1_name) +
        ggtitle(wrapper(paste("Significant kinases in common - ", contrast)))
      save_uka_plot(ukaplot, paste0("./results/common_kinases_", contrast, ".png"))
    }
  }
}

# New helper function to process two UKA data.frames by group (contrast)
overlap_ukas <- function(uka1, uka2, fscorecutoff, data1_name = NULL, data2_name = NULL, prnum = NULL) {
  # Ensure both have Sgroup_contrast
  if (!("Sgroup_contrast" %in% colnames(uka1)) || !("Sgroup_contrast" %in% colnames(uka2))) {
    warning("Sgroup_contrast column not found in one or both UKA files.")
    return(list(commondf = NULL, summary_stats = NULL))
  }
  contrasts <- intersect(unique(as.character(uka1$Sgroup_contrast)), unique(as.character(uka2$Sgroup_contrast)))
  commondf_list <- list()
  summary_stats <- data.frame()
  for (contrast in contrasts) {
    # Parse and filter each group
    parsed1 <- parse_uka(uka1[as.character(uka1$Sgroup_contrast) == contrast, ], fscorecutoff)
    parsed2 <- parse_uka(uka2[as.character(uka2$Sgroup_contrast) == contrast, ], fscorecutoff)
    # Compose comparison string as "prnum - contrast" if prnum is provided
    comparison_str <- if (!is.null(prnum)) paste(prnum, "-", contrast) else contrast
    kinases_info <- get_common_and_unique_kinases(parsed1, parsed2, comparison = comparison_str, data1_name = data1_name, data2_name = data2_name)
    summary_stats <- rbind(summary_stats, kinases_info$summary_row)
    commonids <- kinases_info$commonids
    if (length(commonids) > 0) {
      commondf <- collect_commondf(parsed1, parsed2, commonids, extra_cols = list(comparison = comparison_str))
      commondf_list[[contrast]] <- commondf
    }
  }
  commondf_all <- dplyr::bind_rows(commondf_list)
  list(commondf = commondf_all, summary_stats = summary_stats)
}

process_all_vs_all_uka_summary <- function(uka_files1, uka_files2, fscorecutoff, data1_name, data2_name) {
  all_commondf <- list()
  summary_stats <- data.frame()
  for (i in seq_along(uka_files1)) {
    # Extract prnum from filename (e.g. "180-310_PTK.txt" -> "180-310")
    prnum <- sub("^([^-_]+-[^-_]+)_.*", "\\1", basename(uka_files1[[i]]))
    uka1 <- read_delim(uka_files1[[i]], show_col_types = FALSE) %>% clean_tercen_columns()
    uka2 <- read_delim(uka_files2[[i]], show_col_types = FALSE) %>% clean_tercen_columns()
    result <- overlap_ukas(uka1, uka2, fscorecutoff, data1_name = data1_name, data2_name = data2_name, prnum = prnum)
    if (!is.null(result$commondf) && nrow(result$commondf) > 0) {
      all_commondf[[i]] <- result$commondf
    }
    if (!is.null(result$summary_stats) && nrow(result$summary_stats) > 0) {
      summary_stats <- rbind(summary_stats, result$summary_stats)
    }
  }
  all_commondf_df <- dplyr::bind_rows(all_commondf)
  if (nrow(all_commondf_df) > 0) {
    # Rename columns for CSV output
    colnames(all_commondf_df)[colnames(all_commondf_df) == "data1"] <- data1_name
    colnames(all_commondf_df)[colnames(all_commondf_df) == "data2"] <- data2_name
    all_commondf_df$difference <- all_commondf_df[[data2_name]] - all_commondf_df[[data1_name]]
    ukaplot <- ggplot(all_commondf_df, aes_string(x = data1_name, y = data2_name)) +
      geom_point() +
      facet_wrap(~ scoretype, scales = 'free') +
      ylab(data2_name) +
      xlab(data1_name) +
      ggtitle(wrapper("Significant kinases in common - all comparisons (all files, all contrasts)"))
    save_uka_plot(ukaplot, "./results/common_kinases_all_comparisons_dotplot.png", width = 20, height = 16)
    write_csv(all_commondf_df, "./results/common_kinases_all_comparisons.csv")
    
    make_per_kin_stats(all_commondf_df, data1_name = data1_name, data2_name = data2_name)
  }
  if (nrow(summary_stats) > 0) {
    # No need to rename columns here anymore
    summary_long <- summary_stats %>%
      tidyr::pivot_longer(
        cols = c(overlap_coef, size_ratio, paste0("n_only_", data1_name), paste0("n_only_", data2_name)),
        names_to = "type", values_to = "count"
      )
    # Set type as a factor with desired order
    summary_long$type <- factor(
      summary_long$type,
      levels = c("overlap_coef", "size_ratio", paste0("n_only_", data1_name), paste0("n_only_", data2_name)),
      labels = c("Overlap coefficient", "List size diff (Size Ratio)", paste0("n_only_", data1_name), paste0("n_only_", data2_name))
    )
    boxplot <- ggplot(summary_long, aes(x = "", y = count)) +
      geom_boxplot() +
      geom_jitter(width = 0.1, alpha = 0.5) +
      facet_wrap(~type, scales = 'free_y') +
      ylab("Value") +
      xlab(paste("Kinase set (", data2_name, " vs ", data1_name, ")", sep = "")) +
      ggtitle("Kinase overlap summary across all comparisons") +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.caption = element_text(hjust = 0)
      ) +
      labs(
        caption = "Overlap coef = 1: all smaller list members are in bigger list;\nsize ratio = 1: equal size"
      )
    save_uka_plot(boxplot, "./results/common_kinases_summary_boxplot.png", width = 12, height = 8)
    write_csv(summary_stats, "./results/common_kinases_summary.csv")
  }
}

get_uka_results <- function(uka_files1, uka_files2, output_type, uka_version, data1_name, data2_name, fscorecutoff){
  if (uka_version == "all_vs_all" && output_type == "summary") {
    process_all_vs_all_uka_summary(uka_files1, uka_files2, fscorecutoff, data1_name, data2_name)
  } else {
    for (i in seq_along(uka_files1)){
      if (uka_version == "single_comp"){
        # get comparison
        comparison <- sub(".*\\/(.+).(txt|csv).*", "\\1", uka_files1[[i]])
        print(paste("Working on UKA file:", comparison))
        # get assay type
        assay_type <- ifelse(grepl("_PTK_", uka_files1[[i]]), "PTK", "STK")

        compare_two_uka(ukapath1 = uka_files1[[i]], ukapath2 = uka_files2[[i]],
                        data1_name = data1_name, data2_name = data2_name, fscorecutoff = fscorecutoff,
                        fname = comparison)
      } else if (uka_version == "all_vs_all") {
        if (output_type == "per_comp") {
          process_all_vs_all_uka(uka_files1[[i]], uka_files2[[i]], fscorecutoff, data1_name, data2_name)
        } else if (output_type == "summary") {
          # handled above
        }
      }
    }
  }
}




