# density plot with fraction of middle extreme genes
# in each sub bin of medium expression range
# plus the fraction of bimodal middle extreme genes 
# and fraction of normal middle extreme genes when splitting shape-wise


# import packages ----
library(dplyr)
library(ggplot2)

# parameters ----

# cancer cohorts
cohorts <- readRDS("../data/cohorts.rds")

# shrinking strategy
shrink_strategy <- "sliding_window"

# significance threshold for multiple testing correction
qt <- 0.2

# method for multiple test correction
correction_method <- "BH"

# significance threshold for log-rank test pairwise comparisons
kmpt <- 0.05

# survival analysis
surv_analysis <- "DSS"

# non-bimodal shape of interest (positively_skewed or normal)
shape <- "normal"


# input ----

# shrinked medium bin sizes
shrink_dir <- file.path(
  "../analyses/counts/new/shrink_medium", 
  shrink_strategy, 
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
    sep = "-")
)

# count directory (shape-wise)
count_shape_wise_dir <- file.path(
  "../analyses/counts/new/shape-wise",
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


# output ----

# output directory for plots
plot_dir <- file.path(
  "../plots/kaplan-meier/shrink_medium", 
  shrink_strategy)


# loop ----

# loop over cohorts
density_plot_l <- lapply(seq(cohorts), function(j) {
  
  cat("Process", cohorts[j], "cohort\n") 
  
  # middle extreme bins
  middle_extreme_bins_path <- file.path(
    shrink_dir, 
    cohorts[j], 
    paste0(
      "middle-extreme_", 
      shrink_strategy,
      "-bins_",
      cohorts[j], 
      ".rds")
  )
  
  # count shape-wise
  count_shape_wise_path <- file.path(
    count_shape_wise_dir,
    cohorts[j], 
    surv_analysis, 
    paste0(
      "counts_",
      cohorts[j],
      "_", 
      surv_analysis, 
      ".rds"))
  
  # check if surv-exp file exists for that specific analysis
  if ( file.exists(middle_extreme_bins_path) &
       file.exists(count_shape_wise_path) ) {
    
    # import datasets ----
    
    # survival + expression
    middle_extreme_bins <- readRDS(middle_extreme_bins_path)
    
    # import counts shape-wise
    count_shape_wise <- readRDS(count_shape_wise_path) 
    
    
    # count wrangling ----
    
    # put count dataframe into a list
    count_l <- list(count_shape_wise)
    
    # wrangle each element inside count_l
    count_sub_l <- lapply(seq(count_l), function(i) {
      
      # do not consider to_disentangle hits
      count_df <- count_l[[i]] %>% 
        filter(
          count_low != "to_disentangle", 
          count_medium != "to_disentangle",
          count_high != "to_disentangle")
      
      # transform count columns into numeric
      count_df <- count_df %>% 
        mutate_at(c(
          "count_low",
          "count_medium",
          "count_high"), as.numeric)
      
      # subset non bimodal genes
      count_non_bimodal <- count_df %>% 
        filter(label != "bimodal")
      
      # subset bimodal genes
      count_bimodal <- count_df %>% 
        filter(label == "bimodal")
      
      # return list
      return(
        setNames(
          object = list(count_non_bimodal, count_bimodal), 
          nm = c("non_bimodal", "bimodal")
        )
      )
    })
    
    # assign names to list
    count_sub_l <- setNames(
      object = count_sub_l, 
      nm = c("shape_wise")
    )
    
    # count bimodal shape-wise
    count_bimodal_shape_wise <- count_sub_l$shape_wise$bimodal
    
    # count non bimodal (either normal or positively skewed depending on shape variable)
    count_non_bimodal <- count_sub_l$shape_wise$non_bimodal
    
    # transform some columns into numeric
    middle_extreme_bins[, 3:ncol(middle_extreme_bins)] <- lapply(
      3:ncol(middle_extreme_bins), function(x) as.numeric(middle_extreme_bins[[x]]))
    
    # pick only shape of interest
    middle_extreme_bins <- middle_extreme_bins[middle_extreme_bins$label == shape, ]
    
    # fraction of genes where we see extreme prognosis in each medium sub bin
    fraction_medium_extreme <- sapply(3:ncol(middle_extreme_bins), function(i) {
      
      fraction <- mean(middle_extreme_bins[[i]] %in% c(-1, 1))
      
    })
    
    # transform into dataframe
    fraction_medium_df <- data.frame(
      "bin" = seq(1:length(fraction_medium_extreme)),
      "fraction_bin_medium_extreme" = fraction_medium_extreme)
    
    
    ## shape-wise 
    
    # fraction of bimodal genes in shape-wise
    fraction_bimodal_shape_wise <- mean(
      count_bimodal_shape_wise$count_medium %in% c(-1, 1))
    
    # fraction of non-bimodal genes in shape-wise
    fraction_non_bimodal <- mean(
      count_non_bimodal[count_non_bimodal$label == shape, ]$count_medium %in% c(-1, 1))
    

    # save plot -----
    
    # density plot
    dp <- fraction_medium_df %>% 
      ggplot( aes(x = fraction_bin_medium_extreme)) +
      geom_density(fill = "#69b3a2", color = "#e9ecef", alpha = 0.8) +
      theme_bw() + 
      geom_vline(
        aes(xintercept = fraction_bimodal_shape_wise, color = "bimodal_shape_wise"), 
        linetype = "dashed") +
      geom_vline(
        aes(xintercept = fraction_non_bimodal, color = "non_bimodal"), 
        linetype = "dashed") + 
      scale_color_manual(
        name = "distribution_type", 
        values = c(
          bimodal_shape_wise = "#023e8a", 
          non_bimodal = "#bc4b51")) +
      labs(x = "freq_middle-extreme", y = "Density") + 
      ggtitle(cohorts[j])
    
    # plot file
    plot_file <- file.path(
      plot_dir, 
      paste0(
        "density_fraction_medium-extreme_contiguous_", 
        cohorts[j], 
        ".svg"))
    
    # save
    ggsave(
      filename = plot_file,
      plot = dp,
      device = "svg",
      width = 13.9,
      height = 5.98,
      units = "in"
    )
    
  } else {
    
    cat(
      "No survival-expression file or pvalue file", 
      cohorts[j], " is there\n")
    
    # dataframe with times and counts
    dp <- NULL
    
  } 
  
  return(dp)
  
})


# assign cohorts names to list
density_plot_l <- setNames(
  object = density_plot_l, 
  nm = cohorts)

# cohorts where we do not observe bimodal middle-extreme genes
cohorts_sub <- c(
  "BLCA",
  "CESC", 
  "CHOL",
  "COAD", 
  "DLBC",
  "ESCA", 
  "KICH", 
  "LAML", 
  "LUAD", 
  "PAAD",
  "PCPG",
  "PRAD",
  "READ", 
  "SKCM", 
  "TGCT", 
  "THCA", 
  "THYM", 
  "UCS")

# keep only cohorts where we observe bimodal middle-extreme genes
density_plot_sub_l <- density_plot_l[!names(density_plot_l) %in% cohorts_sub]

# arrange plots together
patchwork::wrap_plots(density_plot_sub_l, ncol = 2, guides = "collect")
