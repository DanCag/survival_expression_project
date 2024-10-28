# +++++++++++++ #
# shrink medium #
# +++++++++++++ #

# split genes in positively skewed and normal distributions
# ensuring the same number of patients as for medium expression group 
# of bimodal genes for each specific cohort. 
# Start from two middle quartiles interval and make sliding windows of size n
# where n is the average medium bin size of bimodal genes in each specific cohort

# the output of this script is needed to produce both figure 3.7 and 3.8
# run the script on workstation

# import packages ----
library(parallel)
library(dplyr)
source("./functions/split_function.R")
source("./functions/km_wrapper.R")
source("./functions/compare_surv-prob.R")
source("./functions/count-prognosis.R")


# parameters ----

# cancer cohorts
cohorts <- readRDS("../data/cohorts.rds")

# significance threshold for multiple testing correction
qt <- 0.2

# method for multiple test correction
correction_method <- "BH"

# significance threshold for log-rank test pairwise comparisons
kmpt <- 0.05

# survival analysis
surv_analysis <- "DSS"

# number of cores for parallelization
n_cores <- 20


# input ----

# directory with survival-expression info
surv_exp_dir <- "../data/survival-tables/process/surv-time_events-censored_no-imp/survival-exp/continuos_exp"

# directory with pvalue of all genes
pvalue_dir <- "../analyses/p-value/new/shape-wise/all-genes"

# file with the median bin size of medium group of bimodal genes in each cohort
median_medium_bin_size_bimodal_path <- "../analyses/patients/new/n_patients_medium-bin_cohorts.tsv"


# output ----

# output directory where I store middle extreme readout for each bin
out_count_dir <- file.path(
  "../analyses/counts/new/shrink_medium/sliding_window", 
  paste0(
    as.character(qt), 
    "-qvalue-threshold"), 
  paste0(
    as.character(kmpt), 
    "-km-pairwise-threshold"))


# data ----

# import median_medium_bin_size_bimodal
mmbsb <- read.delim(
  median_medium_bin_size_bimodal_path, 
  stringsAsFactors = F)


# loop ----

# loop over cohorts
count_summary_cohort_l <- lapply(seq(cohorts), function(j) {
  
  cat("Process", cohorts[j], "cohort\n") 
  
  # surv-exp filepath
  surv_exp_path <- file.path(
    surv_exp_dir, 
    cohorts[j], 
    surv_analysis, 
    paste0(
      "surv-time_exp_events-censored_no-imp_", 
      surv_analysis, 
      "_",
      cohorts[j], 
      ".rds"))
  
  # pvalue path
  pvalue_path <- file.path(
    pvalue_dir, 
    cohorts[j], 
    surv_analysis, 
    paste0(
      "p-value_",
      cohorts[j],
      "_", 
      surv_analysis,
      ".rds"))
  
  # check if surv-exp file exists for that specific analysis
  if ( file.exists(surv_exp_path) & 
       file.exists(pvalue_path) ) {
    
    # import datasets ----
    
    # survival + expression
    surv_exp <- readRDS(surv_exp_path)
    
    # pvalue of all genes
    pvalue <- readRDS(pvalue_path)
    
    
    # survival exp wrangle ----
    
    # make a copy of surv_exp without genes
    # To this df I will add columns with the 
    # expression categories for each gene.
    # This is faster compared with adding the categorical
    # cols to surv_exp which is already a very big table.
    # I need to have survival info and strata info together
    # for later making the survfit object
    surv <- surv_exp[, 1:7]
    
    
    # pvalue wrangling ----
    
    # consider only positively skewed and normal genes with pvalue != NA
    # since NA is assigned when there is no medium group
    genes <- pvalue %>%
      na.omit(pvalue) %>% 
      filter(label %in% c("positive_skewed", "normal")) %>% 
      pull(gene)
    
    # median medium bin size bimodal for cohorts[j] 
    n <- mmbsb %>% 
      filter(cohort == cohorts[j]) %>% 
      pull(median_medium_bin_size)
    
    
    # KM - pvalue ----
    
    # in this first step whre we compute only KM pvalue to later filter for genes
    # with significant pvalues
    
    # loop over genes, fit KM object and extract pvalue of strata comparison
    pval_l <- mclapply(
      genes,
      function(g) {
        
        # cat("Process", g, "\n")
        
        # compute quartiles of gene expression distribution 
        quartiles <- quantile(surv_exp[[g]])
        
        # order gene expression distribution
        distr_ordered <- surv_exp[[g]][order(surv_exp[[g]])]
        
        # subset gene expression distribution keeping only values
        # between 1st quartile and 3rd quartile
        middle_values <- distr_ordered[
          distr_ordered > quartiles[2] & 
            distr_ordered <= quartiles[4]]
        
        # loop over contiguous bins of size n and define the values
        # in the medium bin for each cycle
        medium_l <- lapply(seq(length(middle_values)-(n-1)), function(idx) {
          
          # medium range
          medium <- middle_values[idx:(idx + (n-1))]
          
          # extremes of medium range
          extreme <- medium[c(1, length(medium))]
          
          
          return(extreme)
        })
        
        # loop over medium extremes
        # and compute pvalue for each possible medium group
        pv_v <- sapply(seq(medium_l), function(i) {
          
          # label gene expression
          gene_labels <- ifelse(
            surv_exp[[g]] < medium_l[[i]][1], "low", 
            ifelse(surv_exp[[g]] > medium_l[[i]][2], "high", "medium"))
          
          # pvalue
          pv <- km_pvalue_f(
            gene = g, 
            gene_labels = gene_labels, 
            surv_df = surv)[1]
          
        })
        
        # get the mean pvalue 
        # average_pv <- mean(pv_v)
        min_pv <- min(pv_v)
        
        # result
        result <- c(
          g,
          pvalue$label[
            pvalue$gene == g],
          min_pv)
        
        # return
        return(result)
        
      }, mc.cores = n_cores)
    
    # dataframe with gene, label and pvalue information 
    pval_df <- data.frame(
      "gene" = sapply(pval_l, "[[", 1), 
      "label" = sapply(pval_l, "[[", 2), 
      "min_pvalue" = as.numeric(sapply(pval_l, "[[", 3)))
    
    
    # p-value correction ----
    
    # pvalue correction on testable genes
    # correct shape-wise
    pval_df <- pval_df %>% 
      group_by(label) %>% 
      mutate(
        p_adj = p.adjust(
          min_pvalue, 
          method = correction_method), 
        .after = min_pvalue)
    
    # subset genes that pass multiple test correction
    pval_adj_sig_df <- pval_df[pval_df$p_adj < qt, ]
    
    # check if there is any gene left after mtc
    if ( dim(pval_adj_sig_df)[1] > 0 ) {
      
      # and sort them according to increasing q value
      pval_adj_sig_df <- pval_adj_sig_df[order(pval_adj_sig_df$p_adj), ]
      
      # survival time and probability ----
      
      # loop over these genes and record sub bins where medium 
      # is associated with extreme prognostic behavior
      extreme_medium_bins_l <- mclapply(pval_adj_sig_df$gene, function(cg) {
        
        # cat("Process", cg, "\n")
        
        # compute quartiles of gene expression distribution 
        quartiles <- quantile(surv_exp[[cg]])
        
        # order gene expression distribution
        distr_ordered <- surv_exp[[cg]][order(surv_exp[[cg]])]
        
        # subset gene expression distribution keeping only values
        # between 1st quartile and 3rd quartile
        middle_values <- distr_ordered[
          distr_ordered > quartiles[2] & 
            distr_ordered <= quartiles[4]]
        
        # loop over contiguous bins of size n and define the values
        # in the medium bin for each cycle
        medium_l <- lapply(seq(length(middle_values)-(n-1)), function(idx) {
          
          # medium range
          medium <- middle_values[idx:(idx + (n-1))]
          
          # extremes of medium range
          extreme <- medium[c(1, length(medium))]
          
          return(extreme)
        })
        
        # loop over medium extremes
        # and count if medium is associated 
        # with extreme prognostic behavior in each possible medium group
        counts_l <- lapply(seq(medium_l), function(i) {
          
          # label gene expression
          gene_labels <- factor(
            ifelse(
              surv_exp[[cg]] < medium_l[[i]][1], "low", 
              ifelse(surv_exp[[cg]] > medium_l[[i]][2], "high", "medium")), 
            levels = c("low", "medium", "high"))
          
          # pvalue
          pval <- km_pvalue_f(
            gene = cg, 
            gene_labels = gene_labels, 
            surv_df = surv)[1]
          
          if ( pval < kmpt) {
            
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
              threshold = kmpt)
            
          } else {
            
            # set NAs to counts
            counts_v <- setNames(
              object = rep(NA, 9), 
              nm = c("time_low", "time_medium", "time_high", 
                     "surv_low", "surv_medium", "surv_high", 
                     "count_low", "count_medium", "count_high"))
          }
          
          return(c("bin" = i, counts_v))
          
        })
        
        # dataframe with times and counts
        counts_df <- as.data.frame(do.call(rbind, counts_l))
        
        # add gene name and label to dataframe
        counts_df <- data.frame(
          "gene" = rep(cg, length(medium_l)),
          "label" = rep(pval_adj_sig_df$label[pval_adj_sig_df$gene == cg], length(medium_l)),
          counts_df)
        
        # vector with count_medium information
        extreme_medium <- setNames(
          object = c(
            cg, 
            pval_adj_sig_df$label[pval_adj_sig_df$gene == cg], 
            counts_df$count_medium), 
          nm = c("gene", "label", counts_df$bin))
        
        return(extreme_medium)
        
      }, mc.cores = n_cores) 
      
      # check the number of sub bins obtained for each gene
      n_bins_gene <- sapply(extreme_medium_bins_l, function(x) length(x))
      
      # pick the number of sub bins obtained for most of the genes
      n_bins_max <- names(table(n_bins_gene))[table(n_bins_gene) == max(table(n_bins_gene))]
      
      # keep only genes with the most frequent number of sub bins
      # this is an approximation because we are not considering genes with 
      # different number of sub bins
      elements_to_keep <- sapply(extreme_medium_bins_l, function(x) length(x)) == n_bins_max
      
      # subset list
      extreme_medium_bins_sub_l <- extreme_medium_bins_l[elements_to_keep]
      
      # dataframe
      extreme_medium_bins_df <- as.data.frame(do.call(rbind, extreme_medium_bins_sub_l))
      
      
      # save ----
      
      # output cohort directory
      out_count_cohort_dir <- file.path(
        out_count_dir, 
        cohorts[j])
      
      # make output cohort directory
      dir.create(
        out_count_cohort_dir, 
        showWarnings = F, 
        recursive = T)
      
      # save
      saveRDS(
        extreme_medium_bins_df, 
        file = file.path(
          out_count_cohort_dir, 
          paste0(
            "middle-extreme_sliding_window-bins_",
            cohorts[j],
            ".rds"
          )
        ))
      
    } else {
      
      print("No gene passes the correction\n")
      
      # dataframe with middle-extreme record for all bins
      extreme_medium_bins_df <- NULL
      
    }
    
  } else {
    
    cat(
      "No survival-expression file or pvalue file", 
      cohorts[j], " is there\n")
    
    # dataframe with times and counts
    extreme_medium_bins_df <- NULL
    
  } 
  
  return(extreme_medium_bins_df)
  
})