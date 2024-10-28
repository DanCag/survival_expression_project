# group expression distribution of genes into low, medium 
# and high expression range


# split gene expression of non-bimodal genes into three bins
split_f <- function(gene, split_strategy,  surv_exp_df) {
  
  # cut expression into three bins
  if (split_strategy == "tertiles") {
    
    # "low" == up to 33.3%,
    # "medium" == up to 66.6%
    # and "high" == up to 100%
    
    # bin expression in three groups
    gene_labels <- cut(
      surv_exp_df[[gene]],
      breaks = unique(
        quantile(
          surv_exp_df[[gene]],
          probs = seq.int(0,1,by = 1/3))),
      include.lowest = TRUE)
    
  } else if (split_strategy == "quartiles") {
    
    # "low" == 1st quartile,
    # "medium" == 2nd and 3rd quartiles
    # and "high" == 4th quartile
    
    # bin expression in four groups
    gene_labels <- cut(
      surv_exp_df[[gene]],
      breaks = unique(
        quantile(surv_exp_df[[gene]])),
      include.lowest = TRUE)
  }
  
  # In how many groups the expression values have been cut? (max is 4)
  if ( length(levels(gene_labels)) == 1 ) { # if only 1 group is there
    
    # label the group
    levels(gene_labels) <- "same"
    
  } else if ( length(levels(gene_labels)) == 2 ) { # if 2 groups
    
    # label the groups
    levels(gene_labels) <- c("low", "high")
    
  } else if ( length(levels(gene_labels)) == 3 ) { # if 3 groups
    
    # label the groups
    levels(gene_labels) <- c("low", "medium", "high")
    
  } else { # if 4 groups
    
    # label the groups
    levels(gene_labels) <- c("low", "medium", "medium", "high")
  }
  
  return(gene_labels)
}


# split bimodal gene expression distribution into
# - samples falling into first mode (low)
# - samples falling into interval where minimum between modes is found (medium)
# - samples falling into second mode (high)
split_bimodal3_f <- function(gene, surv_exp_df) {
  
  # spot modes with laplacesdemon package
  modes <- LaplacesDemon::Modes(x = surv_exp_df[[gene]])$modes
  
  # sort modes based on expression value
  # by default values are sorted based on decreasing density of modes
  modes_sorted <- sort(modes)
  
  # compute density function
  d <- density(
    surv_exp_df[[gene]])
  
  # find antimode: minimum value between modes
  antimode <- sapply(seq(length(modes_sorted)-1), function(m) {
    
    # interval between two modes
    interval_bool <- d$x <= modes_sorted[m+1] & d$x >=  modes_sorted[m]
    
    # minimum density in that interval
    min_density <- min(d$y[interval_bool])
    
    # expression value corresponding to minimum density in that interval
    min_exp <- d$x[d$y == min_density]
    
  })
  
  # midpoint low
  midpoint_low <- mean(c(modes_sorted[1], antimode))
  
  # midpoint high
  midpoint_high <- mean(c(modes_sorted[2], antimode))
  
  # defines expression ranges
  gene_labels <- factor(ifelse(
    surv_exp_df[[gene]] < midpoint_low, "low", 
    ifelse(
      surv_exp_df[[gene]] > midpoint_high, 
      "high", "medium")),
    levels = c("low", "medium", "high"))
  
  return(gene_labels)
  
}