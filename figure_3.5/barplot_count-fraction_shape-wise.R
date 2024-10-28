# Plot the count/fraction of genes where one of the three expression group
# is associated with an extreme prognostic behavior (good/bad for survival)


# import packages ----
library(dplyr)
library(ggplot2)
library(svglite)
library(tidyr)


# input variables ----

# survival analysis
surv_analysis <- "DSS"

# cancer cohorts
cohorts <- readRDS("../data/cohorts.rds")

# significance threshold for multiple testing correction
qt <- 0.2

# significance threshold for log-rank test pairwise comparisons
kmpt <- 0.05

# strata names
strata_names <- c("low", "medium", "high")

# directory info
dir_info <- "shape-wise"

# count directory
count_dir <- file.path(
  "../analyses/counts/new",
  dir_info, 
  paste(
    as.character(qt),
    "qvalue", 
    "threshold",
    sep = "-"),
  paste(
    as.character(kmpt),
    "km",
    "pairwise",
    "threshold",
    sep = "-"))

# metrics we want to plot (it could "freq" or "count")
plot_metrics <- "freq"
# plot_metrics <- "count"


# output variables ----

# output directory for plots
output_dir <- file.path(
  "../plots/kaplan-meier/counts/new", 
  dir_info)


# loop ----

# loop across different distribution shapes
count_freq_shape_l <- lapply(c("bimodal", "positive_skewed", "normal"), function(x) {
  
  cat("Process", x, "shape\n")
  
  # loop over cohorts and get a list of count-freq info for each cohort
  count_freq_l <- lapply(seq(length(cohorts)), function(i) {
    
    cat("\n\n", "Process", cohorts[i], "cohort","\n\n",
        "------------------",
        "\n\n")
    
    # import datasets ----
    
    # path of count dataframe
    count_path <- file.path(
      count_dir,
      cohorts[i], 
      surv_analysis, 
      paste0(
        "counts_",
        cohorts[i],
        "_", 
        surv_analysis, 
        ".rds"))
    
    # check if file exists
    if( file.exists(count_path) ) { 
      
      # import counts
      count_df <- readRDS(count_path) 
      
      # wrangling ----
      
      # subset genes with specific shape
      count_shape <- count_df %>% 
        filter(label == x)
      
      # exclude "to_disentangle" hits 
      count_shape_sub <- count_shape %>% 
        filter(
          count_low != "to_disentangle", 
          count_medium != "to_disentangle", 
          count_high != "to_disentangle")
      
      # transform count columns into numeric
      count_shape_sub <- count_shape_sub %>% 
        mutate_at(c("count_low", "count_medium", "count_high"), as.numeric)
      
      # dataframe with cohort, number of genes and number of genes with extreme prognosis
      cohort_info <- data.frame(
        "cohort" = cohorts[i],
        "n_genes" = dim(count_shape_sub)[1],
        "n_genes_extreme_prognosis" = sum(
          table(count_shape_sub$count_low)[names(table(count_shape_sub$count_low)) != "0"],
          table(count_shape_sub$count_medium)[names(table(count_shape_sub$count_medium)) != "0"],
          table(count_shape_sub$count_high)[names(table(count_shape_sub$count_high)) != "0"]))
      
      # count and frequency
      count_freq <- count_shape_sub %>% 
        summarise(
          n_low_bad = sum(count_shape_sub[[paste0("count_", strata_names[1])]] == -1),
          n_low_good = sum(count_shape_sub[[paste0("count_", strata_names[1])]] == 1),
          n_medium_bad = sum(count_shape_sub[[paste0("count_", strata_names[2])]] == -1),
          n_medium_good = sum(count_shape_sub[[paste0("count_", strata_names[2])]] == 1),
          n_high_bad = sum(count_shape_sub[[paste0("count_", strata_names[3])]] == -1),
          n_high_good = sum(count_shape_sub[[paste0("count_", strata_names[3])]] == 1),
          freq_low_bad = n_low_bad/cohort_info$n_genes_extreme_prognosis,
          freq_low_good = n_low_good/cohort_info$n_genes_extreme_prognosis,
          freq_medium_bad = n_medium_bad/cohort_info$n_genes_extreme_prognosis,
          freq_medium_good = n_medium_good/cohort_info$n_genes_extreme_prognosis,
          freq_high_bad = n_high_bad/cohort_info$n_genes_extreme_prognosis,
          freq_high_good = n_high_good/cohort_info$n_genes_extreme_prognosis)
      
      # bind cohort_info and count_freq
      cohort_count_freq <- cbind(cohort_info, count_freq)
      
      
    } else {
      
      print("File does not exist\n\n")
      
      cohort_count_freq <- NULL
      
    }
    
    return(cohort_count_freq)
    
  })
  
  
  # Data wrangling ----
  
  # rbind list of dataframes
  cohort_count_freq_df <- bind_rows(count_freq_l)
  
  # remove cohorts where no expression range is associated with extreme prognosis
  cohort_count_freq_extreme_df <- cohort_count_freq_df[
    cohort_count_freq_df$n_genes_extreme_prognosis != 0, ]
  
  if ( plot_metrics == "count") {
    
    # subset only cohort and count columns and melt the dataframe
    count_melt <- cohort_count_freq_extreme_df %>% 
      select(cohort, 
             n_low_bad,
             n_low_good,
             n_medium_bad, 
             n_medium_good, 
             n_high_bad, 
             n_high_good) %>%
      gather("category", "count", -cohort)
    
    # remove part of string of category col 
    # and transform it into factor
    count_melt$category <- factor(
      sub('.+?_', "", count_melt$category), 
      levels = c("low_bad",
                 "low_good",
                 "medium_bad", 
                 "medium_good",
                 "high_bad", 
                 "high_good"))
    
  } else {
    
    # subset only cohort and frequencies columns and melt the dataframe
    freq_melt <- cohort_count_freq_extreme_df %>% 
      select(cohort, 
             freq_low_bad,
             freq_low_good,
             freq_medium_bad, 
             freq_medium_good, 
             freq_high_bad, 
             freq_high_good) %>%
      gather("category", "freq", -cohort)
    
    
    # remove part of string of category col 
    # and transform it into factor
    freq_melt$category <- factor(
      sub('.+?_', "", freq_melt$category), 
      levels = c("low_bad",
                 "low_good",
                 "medium_bad", 
                 "medium_good",
                 "high_bad", 
                 "high_good"))
  }
  
  
  # Plot ----
  
  # plot dir
  plot_dir <- file.path(
    output_dir, 
    paste(
      as.character(qt), 
      "qvalue", 
      "threshold",
      sep = "-"),
    paste(
      as.character(kmpt), 
      "km", 
      "pairwise", 
      "threshold",
      sep = "-"), 
    x
  )
  
  # create directory(ies) recursively if not there
  dir.create(
    plot_dir, 
    showWarnings = FALSE ,
    recursive = T)
  
  # bar plot
  
  ## stacked barchart
  
  # - count
  if (plot_metrics == "count") {
    
    # barchart
    bar_p <- ggplot(data = count_melt) +
      geom_bar(mapping = aes(x = cohort,
                             y = count, 
                             fill = category),
               stat = "identity",
               position = "stack") +
      labs(title = surv_analysis) +
      theme_bw() + # change legend labels and add no. samples in each cohort
      scale_fill_discrete(
        name = "Expression group",
        labels = c(
          "low-group_bad",
          "low-group_good",
          "medium-group_bad", 
          "medium-group_good", 
          "high-group_bad", 
          "high-group_good")) +
      geom_text(
        data = cohort_count_freq_extreme_df, 
        aes(x = cohort, y = 0, label = n_genes_extreme_prognosis))
    
    # plot file
    plot_file <- file.path(
      plot_dir, 
      paste0(
        "barplot_count_", 
        x, 
        "_BH-significant",
        paste0("-", as.character(qt)),
        "_",
        surv_analysis, 
        ".svg"))
    
    # save plot
    # plotting like this, 
    # it does not plot the number of genes in each cohort
    ggsave( 
      filename = plot_file,
      plot = bar_p, 
      device = "svg", 
      width = 13.9, 
      height = 5.98, 
      units = "in"
    )
    
  } else { # -freq
    
    # barchart
    bar_p <- ggplot(data = freq_melt) +
      geom_bar(mapping = aes(x = cohort,
                             y = freq, 
                             fill = category),
               stat = "identity",
               position = "stack") +
      labs(title = surv_analysis) +
      theme_bw() +
      scale_fill_discrete(
        name = "Expression group",
        labels = c(
          "low-group_bad",
          "low-group_good",
          "medium-group_bad", 
          "medium-group_good", 
          "high-group_bad", 
          "high-group_good"))+ 
      geom_text(
        data = cohort_count_freq_extreme_df,
        aes(x = cohort, y = 0, label = n_genes_extreme_prognosis))
    
    # plot file
    plot_file <- file.path(
      plot_dir, 
      paste0(
        "barplot_freq_", 
        x, 
        "_BH-significant",
        paste0("-", as.character(qt)),
        "_",
        surv_analysis, 
        ".svg")
    )
    
    # save
    ggsave( 
      filename = plot_file,
      plot = bar_p, 
      device = "svg", 
      width = 13.9, 
      height = 5.98, 
      units = "in"
    )
    
    # plot where good/bad prognosis are merged into extreme prognosis
    
    # make a copy of freq_melt
    freq_melt2 <- freq_melt
    
    # make a macro category column
    freq_melt2$category_macro <- factor(
      sub("_.*","",freq_melt2$category), 
      levels = c("low", "medium", "high"))
    
    # barplot where x and y axis are flipped
    bar_p_macro <- ggplot(data = freq_melt2) +
      geom_bar(mapping = aes(x = forcats::fct_rev(cohort), y = freq, fill = category_macro), stat = "identity",
               position = "stack") +
      scale_fill_manual(values = c("#0a9396ff", "#ee9b00ff", "#bb3e03ff")) + 
      labs(title = surv_analysis, y = "frequency", x = "cohort") +
      theme_bw() +
      coord_flip()+
      # scale_fill_discrete(
      #   name = "Expression group")+ 
      geom_text(
        data = cohort_count_freq_extreme_df,
        aes(x = cohort, y = 0, label = n_genes_extreme_prognosis))
    
    # plot file
    plot_macro_file <- file.path(
      plot_dir, 
      paste0(
        "barplot_freq-macro_", 
        x, 
        "_BH-significant",
        paste0("-", as.character(qt)),
        "_",
        surv_analysis, 
        ".svg")
    )
    
    # save
    ggsave( 
      filename = plot_macro_file,
      plot = bar_p_macro, 
      device = "svg")
    
  }
  
})