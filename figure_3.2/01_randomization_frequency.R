# +++++++++++++ #
# Randomization #
# +++++++++++++ #

# Shuffle gene expression keeping survival information stable, 
# group by low, medium and high expression range and check if 
# you get genes significantly associated with survival by chance
# Perform N iterations on each cohort. For each iteration register
# the number of genes in low, medium and high significantly associated with survival

# Run this script on workstation
# Runtime on workstation:
# 100 iterations, all cohort, 20 cores, mclapply on pval_l and count_l loops 
# total time --> 18 hours


# start time
start_time <- Sys.time()

# import packages ----
library(parallel)
library(dplyr)
source("./functions/split_function.R")
source("./functions/km_wrapper.R")
source("./functions/compare_surv-prob.R")
source("./functions/count-prognosis.R")


# input variables ----

# survival analysis
surv_analysis <- "DSS"

# cancer cohorts
cohorts <- readRDS("../data/cohorts.rds")

# split by
split_by <- "tertiles"

# significance threshold for multiple testing correction
qvalue_threshold <- 0.2

# method for multiple test correction
method <- "BH"

# number of iterations
n_iterations <- 100

# seeds for iteration
set.seed(123456)

# get a different seed for each iteration to be sure the shuffling changes
# and that you can track back the exact shuffling
seeds <- sample.int(100, n_iterations) 

# number of cores for parallelization
n_cores <- 20


## input

# directory with survival-expression info
surv_exp_dir <- "../data/survival-tables/process/surv-time_events-censored_no-imp/survival-exp/continuos_exp"


## output

# output directory
out_dir <- file.path(
  "../analyses/randomization/count", 
  paste0(
    as.character(n_iterations),
    "_iterations"))

# create directory(ies) recursively if not there
dir.create(
  out_dir,
  showWarnings = FALSE ,
  recursive = T)

# output file
output_file <- file.path(
  out_dir, 
  paste0(
    "count_", 
    as.character(n_iterations),
    "-iterations.rds")
)


# loop ----


# loop over iterations
count_cohort_iter_l <- lapply(seq(cohorts), function(j) {
  
  cat("Process", cohorts[j], "cohort\n") 
  
  # surv-exp filepath
  surv_exp_path <- file.path(
    surv_exp_dir, 
    cohorts[j], 
    surv_analysis, 
    paste0("surv-time_exp_events-censored_no-imp_", 
           surv_analysis, 
           "_",
           cohorts[j], 
           ".rds")
  )
  
  # check if surv-exp file exists for that specific analysis
  if ( file.exists(surv_exp_path) ) {
    
    # import datasets ----
    
    # survival + expression
    surv_exp <- readRDS(surv_exp_path)
    
    # make a copy of surv_exp without genes
    # To this df I will add columns with the 
    # expression categories for each gene.
    # This is faster compared with adding the categorical
    # cols to surv_exp which is already a very big table.
    # I need to have survival info and strata info together
    # for later making the survfit object
    surv <- surv_exp[, 1:7]
    
    
    # Data wrangling ----
    
    # subset cols with gene expression
    cols_exp <- colnames(surv_exp)[
      grep("ENS", colnames(surv_exp))]
    
    # remove genes with all 0s 
    cols_exp_selected <- cols_exp[
      apply(surv_exp[, cols_exp], 
            2, 
            function(x) sum(x) != 0)]
    
    # randomization ----
    
    # loop over n_iterations
    iter_l <- lapply(seq(n_iterations), function(it){
      
      cat("Iteration #", it, "\n\n")
      
      # set seed
      set.seed(seeds[it])
      
      # random permutation of rows
      surv_exp_shuffled <- surv_exp[sample(nrow(surv_exp)), ]
      
      
      # loop ----
      
      # loop over genes, fit KM object and extract pvalue of strata comparison
      pval_l <- mclapply(
        cols_exp_selected,
        function(g) {
          
          # cat("Process", g, "\n")
          
          # gene labels
          gene_labels <- split_f(
            gene = g,
            split_strategy = split_by,
            surv_exp_df = surv_exp_shuffled)
          
          # Is "medium" group there?
          if ( "medium" %in% gene_labels &
               min(table(gene_labels)) >= 3 ) {
            
            # pvalue
            pv <- km_pvalue_f(
              gene = g, 
              gene_labels = gene_labels, 
              surv_df = surv)[1]
            
          } else {
            
            # check if previous if failed because "medium" group is not there
            # or because there were < 3 samples in on either "low", "medium" or
            # "high" group
            
            # medium group is not there 
            if ( !"medium" %in% gene_labels ) {
              
              pv <- NA
              
            } else { # not enough samples in at least one of the three groups
              pv <- "not_enough_statistical_power"
            }
            
          }
          
          return(pv)
        }, mc.cores = n_cores)
      
      
      # multiple testing correction ----
      
      # unlist pval_l and assign genes' names
      pval <- setNames(
        object = unlist(pval_l), 
        nm = cols_exp_selected)
      
      # remove NAs
      pval_no_NA <- na.omit(pval)
      
      # adjust pvalue
      p_adj <- p.adjust(pval_no_NA, method = method)
      
      # subset genes that pass multiple test correction
      p_adj_sig <- p_adj[p_adj < qvalue_threshold]
      
      
      # count ----
      
      # check if there is any gene left after mtc
      if ( length(p_adj_sig) > 0 ) {
        
        # sort it according to increasing q value
        p_adj_sig <- p_adj_sig[order(p_adj_sig)]
        
        # loop over these genes and count genes where either 
        # low, medium or high is associated with extreme prognostic behavior
        counts_l <- mclapply(names(p_adj_sig), function(cg) {
          
          # gene labels
          gene_labels <- split_f(
            gene = cg,
            split_strategy = split_by,
            surv_exp_df = surv_exp_shuffled)
          
          # pairwise pvalues
          pairwise_pvals <- km_pvalue_f(
            gene = cg, 
            gene_labels = gene_labels, 
            surv_df = surv)[2:4]
          
          # fit dataframe ordered
          fit_df <- km_fit_f(
            gene = cg,
            gene_labels = gene_labels,
            surv_df = surv)
          
          # compute survival probability
          # year: start from the max survival available 
          # and in case you don't find any time suitable in three strata
          # decrease by 1 
          time_sp_l <- surv_prob_comp_3_groups(
            fit_df = fit_df,
            cg = paste0(cg, "_cat"),
            year = ceiling(max(fit_df$time)),
            strata_names = sub(".*=", "", levels(fit_df$strata)))
          
          # compare survival probabilities and count prognostic behaviors
          # threshold is the significance threshold for pairwise KM comparisons
          counts_v <- counts_3_groups(
            time_sp_l,
            pairwise_info = pairwise_pvals,
            strata_names = sub(".*=", "", levels(fit_df$strata)), 
            threshold = 0.05)
          
          return(counts_v)
          
        }, mc.cores = n_cores)
        
        # dataframe with times and counts
        counts_df <- as.data.frame(do.call(rbind, counts_l))
        
        # add gene name and label to dataframe
        counts_df <- data.frame(
          "gene" = names(p_adj_sig),
          counts_df)
        
        # remove genes with NAs
        # NAs might come from the fact that the surv_prob_comp function
        # does not manage compare survival probs at too much different times
        counts_df <- counts_df[rowSums(is.na(counts_df)) != 9, ]
        
        # number of genes where either low, medium or high show
        # extreme behavior
        low_bad_sig <- sum(counts_df$count_low == -1)
        freq_low_bad_sig <- sum(counts_df$count_low == -1)/length(cols_exp_selected)
        low_good_sig <- sum(counts_df$count_low == 1)
        freq_low_good_sig <- sum(counts_df$count_low == 1)/length(cols_exp_selected)
        medium_bad_sig <- sum(counts_df$count_medium == -1)
        freq_medium_bad_sig <- sum(counts_df$count_medium == -1)/length(cols_exp_selected)
        medium_good_sig <- sum(counts_df$count_medium == 1)
        freq_medium_good_sig <- sum(counts_df$count_medium == 1)/length(cols_exp_selected)
        high_bad_sig <- sum(counts_df$count_high == -1)
        freq_high_bad_sig <- sum(counts_df$count_high == -1)/length(cols_exp_selected)
        high_good_sig <- sum(counts_df$count_high == 1)
        freq_high_good_sig <- sum(counts_df$count_high == 1)/length(cols_exp_selected)
        
        # vector with results
        res <- setNames(
          object = c(
            cohorts[j],
            it, 
            length(cols_exp_selected),
            low_bad_sig,
            freq_low_bad_sig,
            low_good_sig, 
            freq_low_good_sig,
            medium_bad_sig, 
            freq_medium_bad_sig,
            medium_good_sig, 
            freq_medium_good_sig,
            high_bad_sig, 
            freq_high_bad_sig,
            high_good_sig, 
            freq_high_good_sig),
          nm = c(
            "cohort", 
            "iteration",
            "n_genes_tested",
            "n_low_bad", 
            "freq_low_bad",
            "n_low_good", 
            "freq_low_good",
            "n_medium_bad",
            "freq_medium_bad",
            "n_medium_good",
            "freq_medium_good",
            "n_high_bad",
            "freq_high_bad",
            "n_high_good",
            "freq_high_good"
          ))
        
      } else {
        
        print("No gene passes the correction\n")
        
        # vector with results
        res <- setNames(
          object = c(
            cohorts[j],
            it,
            length(cols_exp_selected), 
            rep(0, 12)), 
          nm = c(
            "cohort",
            "iteration",
            "n_genes_tested",
            "n_low_bad", 
            "freq_low_bad",
            "n_low_good", 
            "freq_low_good",
            "n_medium_bad",
            "freq_medium_bad",
            "n_medium_good",
            "freq_medium_good",
            "n_high_bad",
            "freq_high_bad",
            "n_high_good",
            "freq_high_good"
          ))
        
      }
      
    })
    
    # pool iteration results
    iter_df <- as.data.frame(do.call(rbind, iter_l))
    
    
  } else {
    
    cat(
      "No survival-expression file or labeled genes for", 
      cohorts[j], " is there\n")
    
    # iteration results
    iter_df <- NULL
    
  }
  
  
  # save ----
  
  # directory where to store iteration result
  cohort_dir <- file.path(
    out_dir, 
    cohorts[j])
  
  # create directory(ies) recursively if not there
  dir.create(
    cohort_dir,
    showWarnings = FALSE ,
    recursive = T)
  
  # save res 
  saveRDS(
    object = iter_df, 
    file = file.path(
      cohort_dir, 
      paste0(
        cohorts[j], 
        "_",
        as.character(n_iterations),
        "-iterations.rds")
    )
  )
  
  return(iter_df)
  
})


# pool results
count_cohort_iter_df <- as.data.frame(
  do.call(rbind, count_cohort_iter_l))

# transform count column into numeric
count_cohort_iter_df <- count_cohort_iter_df %>% 
  mutate_at(
    c(
      "iteration",
      "n_genes_tested",
      "n_low_bad", 
      "freq_low_bad",
      "n_low_good", 
      "freq_low_good",
      "n_medium_bad",
      "freq_medium_bad",
      "n_medium_good",
      "freq_medium_good",
      "n_high_bad",
      "freq_high_bad",
      "n_high_good",
      "freq_high_good"
    ), 
    as.numeric)


# save count_iter_cohort_df
saveRDS(
  count_cohort_iter_df, 
  file = output_file
)

# end time
end_time <- Sys.time()