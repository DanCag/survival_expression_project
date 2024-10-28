# +++++++++++++ #
# Count barplot #
# +++++++++++++ #

# Plot the count/fraction of genes where one of the three expression group
# is associated with an extreme prognostic behavior (good/bad for survival)
# in cancer cohorts 


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

# strategy used to split gene expression in ranges
split_by <- "tertiles"

# significance threshold for multiple testing correction
qvalue_threshold <- 0.2

# significance threshold for log-rank test pairwise comparisons
km_pairwise_threshold <- 0.05

# strata names change depending if 
# we use bimodal or quantile splitting strategy
strata_names <- c("low", "medium", "high")

# count directory
count_dir <- file.path(
  "../analyses/counts/new",
  split_by,
  paste(
    as.character(qvalue_threshold),
    "qvalue", 
    "threshold",
    sep = "-"),
  paste(
    as.character(km_pairwise_threshold),
    "km",
    "pairwise",
    "threshold",
    sep = "-"))

# metrics we want to plot (it can be "fraction" or "count")
# plot_metrics <- "fraction"
plot_metrics <- "count"


# output variables ----

# plot dir
plot_dir <- file.path(
  "../plots/kaplan-meier/counts/new", 
  split_by,
  paste(
    as.character(qvalue_threshold), 
    "qvalue", 
    "threshold",
    sep = "-"),
  paste(
    as.character(km_pairwise_threshold), 
    "km", 
    "pairwise", 
    "threshold",
    sep = "-")
)

# create directory(ies) recursively if not there
dir.create(
  plot_dir, 
  showWarnings = FALSE ,  
  recursive = T)


# loop ----

# loop over cohorts and get a list with count/fraction of genes with extreme prog
count_fraction_l <- lapply(seq(length(cohorts)), function(i) {
  
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
      ".rds")
  )
  
  # check if file exists
  if( file.exists(count_path) ) { 
    
    # counts info
    count_df <- readRDS(count_path)
    
    # Data wrangling ----  
    
    # exclude "to_disentangle" hits 
    count_df_sub <- count_df %>% 
      filter(
        count_low != "to_disentangle", 
        count_medium != "to_disentangle", 
        count_high != "to_disentangle")
    
    # transform count columns into numeric
    count_df_sub <- count_df_sub %>% 
      mutate_at(c("count_low", "count_medium", "count_high"), as.numeric)
    
    # discard genes where there is no effect on survival: 
    # count_low & count_medium & count_high == 0
    genes_to_keep <- rowSums(abs(count_df_sub[, -c(1:7)])) != 0
    count_effect <- count_df_sub[genes_to_keep, ]
    
    # dataframe 
    cohort_info <- data.frame(
      "cohort" = cohorts[i],
      "n_genes_where_we_could_compare_surv" = dim(count_df_sub)[1],
      "n_genes_extreme_prognosis" = sum(
        table(count_effect$count_low)[names(table(count_effect$count_low)) != 0],
        table(count_effect$count_medium)[names(table(count_effect$count_medium)) != 0],
        table(count_effect$count_high)[names(table(count_effect$count_high)) != 0]), 
      "n_genes_extreme_prognosis_unique" = dim(count_effect)[1])
    
    # count and fraction
    # fraction is computed on the total number of genes with extreme prognosis
    # counting the same gene twice in case 2 expression ranges show extreme prog
    count_fraction <- count_effect %>% 
      summarise(
        n_low_bad = sum(count_effect[[paste0("count_", strata_names[1])]] == -1),
        n_low_good = sum(count_effect[[paste0("count_", strata_names[1])]] == 1),
        n_medium_bad = sum(count_effect[[paste0("count_", strata_names[2])]] == -1),
        n_medium_good = sum(count_effect[[paste0("count_", strata_names[2])]] == 1),
        n_high_bad = sum(count_effect[[paste0("count_", strata_names[3])]] == -1),
        n_high_good = sum(count_effect[[paste0("count_", strata_names[3])]] == 1),
        fraction_low_bad = n_low_bad/cohort_info$n_genes_extreme_prognosis,
        fraction_low_good = n_low_good/cohort_info$n_genes_extreme_prognosis,
        fraction_medium_bad = n_medium_bad/cohort_info$n_genes_extreme_prognosis,
        fraction_medium_good = n_medium_good/cohort_info$n_genes_extreme_prognosis,
        fraction_high_bad = n_high_bad/cohort_info$n_genes_extreme_prognosis,
        fraction_high_good = n_high_good/cohort_info$n_genes_extreme_prognosis)
    
    # bind cohort_info and count_fraction
    cohort_count_fraction <- cbind(cohort_info, count_fraction)
    
  } else {
    
    print("File does not exist\n\n")
    
    # bind cohort_info and count_fraction
    cohort_count_fraction <- NULL
    
  }
  
  return(cohort_count_fraction)
  
})


# Data wrangling ----

# rbind list of dataframes
cohort_count_fraction_df <- bind_rows(count_fraction_l)

# save table with this info
write.table(
  cohort_count_fraction_df[, c(1:10)], 
  file = "../analyses/counts/new/tertiles/0.2-qvalue-threshold/0.05-km-pairwise-threshold/info_all-cohorts.tsv", 
  sep = "\t", 
  row.names = F, 
  col.names = T
)


# remove cohorts where no expression range is associated with extreme prognosis
cohort_count_fraction_extreme_df <- cohort_count_fraction_df[
  cohort_count_fraction_df$n_genes_extreme_prognosis_unique != 0, ]

if ( plot_metrics == "count") {
  
  # subset only cohort and count columns and melt the dataframe
  count_melt <- cohort_count_fraction_extreme_df %>% 
    select(
      cohort, 
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
  
  # subset only cohort and fractionuencies columns and melt the dataframe
  fraction_melt <- cohort_count_fraction_extreme_df %>% 
    select(
      cohort, 
      fraction_low_bad,
      fraction_low_good,
      fraction_medium_bad, 
      fraction_medium_good, 
      fraction_high_bad, 
      fraction_high_good) %>%
    gather("category", "fraction", -cohort)
  
  # remove part of string of category col 
  # and transform it into factor
  fraction_melt$category <- factor(
    sub('.+?_', "", fraction_melt$category), 
    levels = c("low_bad",
               "low_good",
               "medium_bad", 
               "medium_good",
               "high_bad", 
               "high_good"))

}


# Plot ----
  
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
      data = cohort_count_fraction_extreme_df, 
      aes(x = cohort, y = 0, label = n_genes_extreme_prognosis))
  
  # plot file
  plot_file <- file.path(
    plot_dir, 
    paste0(
      "barplot_count_BH-significant",
      paste0("-", as.character(qvalue_threshold)),
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
  
} else { # -fraction
  
  # barchart
  bar_p <- ggplot(data = fraction_melt) +
    geom_bar(mapping = aes(x = cohort,
                           y = fraction, 
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
      data = cohort_count_fraction_extreme_df,
      aes(x = cohort, y = 0, label = n_genes_extreme_prognosis))
  
  # plot file
  plot_file <- file.path(
    plot_dir, 
    paste0(
      "barplot_fraction_BH-significant",
      paste0("-", as.character(qvalue_threshold)),
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
  
  # make a copy of fraction_melt
  fraction_melt2 <- fraction_melt
  
  # make a macro category column
  fraction_melt2$category_macro <- factor(
    sub("_.*","",fraction_melt2$category), 
    levels = c("low", "medium", "high"))
  
  # barplot where x and y axis are flipped
  bar_p_macro <- ggplot(data = fraction_melt2) +
    geom_bar(mapping = aes(x = forcats::fct_rev(cohort), y = fraction, fill = category_macro), stat = "identity",
             position = "stack") +
    scale_fill_manual(values = c("#0a9396ff", "#ee9b00ff", "#bb3e03ff")) + 
    labs(title = surv_analysis, y = "fraction", x = "cohort") +
    theme_bw() +
    coord_flip() +
    # scale_fill_discrete(
    #   name = "Expression group")+ 
    geom_text(
      data = cohort_count_fraction_extreme_df,
      aes(x = cohort, y = 0, label = n_genes_extreme_prognosis))
  
  # plot file
  plot_macro_file <- file.path(
    plot_dir, 
    paste0(
      "barplot_fraction-macro_BH-significant",
      paste0("-", as.character(qvalue_threshold)),
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