# classify distributions
# group by low, medium and high expression range and record 
# p-value of genes significantly associated with survival 
# for each gene that passes multiple testing correction, 
# compare survival time and define which expression range is associated
# with extreme prognosis 

# fig 3.1b obtained splitting by tertiles irrespectively of shape
# fig 3.5 obtained splitting by quartiles (normal-like and positively-skewed disrtibutions)
# and splitting bimodal distributions considering as "medium" the area surrounding
# local minimum between the two modes


# Run this script on workstation

# start time
start_time <- Sys.time()

# import packages ----
library(parallel) # package for parallelization
library(dplyr) # package for data wrangling
source("./functions/classification_f.R")
source("./functions/split_function.R")
source("./functions/km_wrapper.R")
source("./functions/compare_surv-prob.R")
source("./functions/count-prognosis.R")


# parameters ----

# cancer cohorts
cohorts <- readRDS("../data/cohorts.rds")

# split by 
# strategy to split non bimodal distribution
# it can be either "tertiles" or "quartiles"
split_by <- "tertiles"
# split_by <- "quartiles"

# variable to decide whether to split based on shape or not
# split_based_on_shape <- T
split_based_on_shape <- F


# define directory_info depending on split_based_on_shape
if ( split_based_on_shape ) {
  
  directory_info <- "shape-wise"
  
  } else {
    
    if (split_by == "tertiles") {
      
      directory_info <- "tertiles"
      
    } else {
      
      directory_info <- "quartiles"
      
    }
    
  }

# minimum number of patients per group (low/medium/high)
min_n_patients <- 3

# bimodal threshold
bt <- 0.05

# normality threshold
nt <- 0.05

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


# output ----

# output directory where I store classification info
out_cl_dir <- "../analyses/kaplan-meier/classification/TCGA"

# output directory where I store untestable genes
# - gene with no medium group
# - gene with < N patients in one of the three expression range
out_untestable_dir <- file.path(
  "../analyses/genes_list/new", 
  directory_info, 
  "untestable")

# output directory where I store p-value
out_pval_dir <- file.path(
  "../analyses/p-value/new", 
  directory_info, 
  "all-genes")

# output directory where I store genes that don't pass mtc
out_non_sig_dir <- file.path(
  "../analyses/genes_list/new", 
  directory_info, 
  "non-significant")

# output directory where I store p-value and q-value of genes that pass mtc
out_sig_dir <- file.path(
  "../analyses/p-value/new", 
  directory_info, 
  "significant-genes")

# output directory where I store time, survival info
out_count_dir <- file.path(
  "../analyses/counts/new", 
  directory_info, 
  paste0(
    as.character(qt), 
    "-qvalue-threshold"), 
  paste0(
    as.character(kmpt), 
    "-km-pairwise-threshold"))


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
    
    
    # classification ----
    
    # loop over genes and classify their distributions
    classification_l <- mclapply(cols_exp_selected, function(g) {
      
      # cat("Process", g, "\n")
      
      # classification
      cl <- classification_f(gene = g, exp = surv_exp)
      
    }, mc.cores = n_cores)
    
    # bind rows
    classification_df <- as.data.frame(do.call(rbind, classification_l))
    
    # directory to store classification info
    out_cl_cohort_dir <- file.path(
      out_cl_dir,
      cohorts[j])

    # make directory to store classification info
    dir.create(
      out_cl_cohort_dir,
      showWarnings = F,
      recursive = T)

    # save information
    saveRDS(
    classification_df,
    file = file.path(
      out_cl_cohort_dir,
      paste0(
        "genes_classification_",
        cohorts[j],
        ".rds"
      )))

    # subset of genes for which shape is identified
    cols_exp_labeled <- classification_df$gene[classification_df$label != "unlabeled"]
    
    
    # KM - pvalue ----
    
    # loop over genes, fit KM object and extract pvalue of strata comparison
    pval_l <- mclapply(
      cols_exp_labeled,
      function(g) {
        
        # cat("Process", g, "\n")
        
        if ( classification_df$label[
          classification_df$gene == g] %in% c("normal", "positive_skewed") ) {
          
          # gene labels
          gene_labels <- split_f(
            gene = g,
            split_strategy = split_by,
            surv_exp_df = surv_exp)
          
          if ( "medium" %in% gene_labels &
               min(table(gene_labels)) >= min_n_patients ) {
            
            # pvalue
            pv <- km_pvalue_f(
              gene = g, 
              gene_labels = gene_labels, 
              surv_df = surv)[1]
            
          } else {
            
            # medium group is not there 
            if ( !"medium" %in% gene_labels ) {
              
              pv <- NA
              
            } else { # not enough samples in at least one of the three groups
              
              pv <- "not_enough_statistical_power"
            }
            
          }
          
        } else { # we deal with bimodal genes
          
          # if we split bimodal genes differently considering antimode and midpoints
          if ( split_based_on_shape) {
            
            # bin expression in three groups
            gene_labels <- split_bimodal3_f(
              gene = g, 
              surv_exp_df = surv_exp)
            
          } else {
            
            # split in quartiles even if distribution is bimodal
            gene_labels <- split_f(
              gene = g,
              split_strategy = split_by,
              surv_exp_df = surv_exp)
          }

          
          if ( "medium" %in% gene_labels &
               min(table(gene_labels)) >= min_n_patients ) {
            
            # pvalue
            pv <- km_pvalue_f(
              gene = g, 
              gene_labels = gene_labels, 
              surv_df = surv)[1]
            
          } else {
            
            # medium group is not there 
            if ( !"medium" %in% gene_labels ) {
              
              pv <- NA
              
            } else { # not enough samples in at least one of the three groups
              
              pv <- "not_enough_statistical_power"
            }
          }
          
        }
        
        
        # result
        result <- c(
          g,
          classification_df$label[
            classification_df$gene == g], 
          pv)
        
        # return
        return(result)
        
      }, mc.cores = n_cores)
    
    
    # save p-value ----
    
    # dataframe with gene, label and pvalue information 
    pval_df <- data.frame(
      "gene" = sapply(pval_l, "[[", 1), 
      "label" = sapply(pval_l, "[[", 2), 
      "pvalue" = sapply(pval_l, "[[", 3))
    
    # save pval
    out_pval_cohort_dir <- file.path(
      out_pval_dir, 
      cohorts[j], 
      surv_analysis)
    
    # create directory(ies) recursively if not there
    dir.create(
      out_pval_cohort_dir,
      showWarnings = FALSE ,
      recursive = T)
    
    # save p-value
    saveRDS(
      object = pval_df, 
      file = file.path(
        out_pval_cohort_dir, 
        paste0(
          "p-value_", 
          cohorts[j], 
          "_",
          surv_analysis, 
          ".rds")
      )
    )
    
    
    # save untestable genes ----
    
    # save untestable genes
    # - gene with no medium group
    # - gene with < N patients in one of the three expression range
    if ( any(is.na(pval_df$pvalue)) | 
         any(pval_df$pvalue == "not_enough_statistical_power") ) {
      
      # untestable genes
      out_untestable_cohort_dir <- file.path(
        out_untestable_dir,
        cohorts[j])
      
      # create directory(ies) recursively if not there
      dir.create(
        out_untestable_cohort_dir,
        showWarnings = FALSE ,
        recursive = T)
      
      # corresponding to genes where pvalue is NA
      write.table(
        x = pval_df$gene[
          (is.na(pval_df$pvalue)) | 
            (pval_df$pvalue == "not_enough_statistical_power")],
        file = file.path(
          out_untestable_cohort_dir,
          paste0(
            "genes_untestable_",
            cohorts[j],
            ".txt")),
        quote = F,
        col.names = F,
        row.names = F)
      
    }
    
    
    # multiple testing correction ----
    
    # remove untestable genes (NA or not_enough_statistical_power)
    pval_testable_df <- pval_df[
      !(is.na(pval_df$pvalue)) & 
        (pval_df$pvalue != "not_enough_statistical_power"), ]
    
    # pvalue correction on testable genes
    # correct shape-wise
    pval_testable_df <- pval_testable_df %>% 
      group_by(label) %>% 
      mutate(
        p_adj = p.adjust(
          pvalue, 
          method = correction_method), 
        .after = pvalue)

    # subset genes that pass multiple test correction
    pval_adj_sig_df <- pval_testable_df[pval_testable_df$p_adj < qt, ]
    
    
    ## save genes that did not pass multiple test correction
    
    # not significant pathway
    out_non_sig_cohort_dir <- file.path(
     out_non_sig_dir,
      paste0(
        as.character(qt),
        "-qvalue-threshold"),
      cohorts[j],
      surv_analysis)
    
    # create directory(ies) recursively if not there
    dir.create(
      out_non_sig_cohort_dir,
      showWarnings = FALSE ,
      recursive = T)
    
    # save directory
    write.table(
      x = pval_testable_df$gene[pval_testable_df$p_adj >= qt],
      file = file.path(
        out_non_sig_cohort_dir,
        paste0(
          "genes_non-significant-after-mtc_", 
          cohorts[j], 
          "_",
          surv_analysis, 
          ".txt")),
      quote = F,
      col.names = F,
      row.names = F
    )
    
    
    # check if there is any gene left after mtc
    if ( dim(pval_adj_sig_df)[1] > 0 ) {
      
      # and sort them according to increasing q value
      pval_adj_sig_df <- pval_adj_sig_df[order(pval_adj_sig_df$p_adj), ]
      
      # output dir
      out_sig_cohort_dir <- file.path(
        out_sig_dir,
        cohorts[j],
        surv_analysis)
      
      # create directory(ies) recursively if not there
      dir.create(
        out_sig_cohort_dir,
        showWarnings = FALSE ,
        recursive = T)
      
      ## save genes that pass multiple test correction
      saveRDS(
        object = pval_adj_sig_df,
        file = file.path(
          out_sig_cohort_dir,
          paste0(
            "p-value_significant-after-mtc_", 
            cohorts[j], 
            "_",
            surv_analysis, 
            ".rds")))
      
      
      # survival time and probability ----
      
      # loop over these genes and count genes where either 
      # low, medium or high is associated with extreme prognostic behavior
      counts_l <- mclapply(pval_adj_sig_df$gene, function(cg) {
        
        if ( pval_adj_sig_df$label[
          pval_adj_sig_df$gene == cg] %in% c("normal", "positive_skewed") ) {
          
          # gene labels
          gene_labels <- split_f(
            gene = cg,
            split_strategy = split_by,
            surv_exp_df = surv_exp)
          
        } else { # we deal with bimodal genes
          
          # if we split bimodal differently considering antimode and midpoints
          if ( split_based_on_shape) {
            
            # bin expression in three groups
            gene_labels <- split_bimodal3_f(
              gene = cg, 
              surv_exp_df = surv_exp)
            
          } else {
            
            # split in quartiles even if distribution is bimodal
            gene_labels <- split_f(
              gene = cg,
              split_strategy = split_by,
              surv_exp_df = surv_exp)
          }
          
        }
        
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
        "gene" = pval_adj_sig_df$gene,
        "label" = pval_adj_sig_df$label,
        counts_df)
      
      # remove genes with NAs
      # NAs might come from the fact that the surv_prob_comp function
      # does not manage compare survival probs at too much different times
      counts_df <- counts_df[rowSums(is.na(counts_df)) != 9, ]
      
      
      # save ----
      
      if ( nrow(counts_df) == 0 ) {
        
        cat("All NAs\n")
        
      } else {
        
        # directory where I save KM results of cohorts[j]
        out_count_cohort_dir <- file.path(
          out_count_dir, 
          cohorts[j], 
          surv_analysis)
        
        # create directory(ies) recursively if not there
        dir.create(
          out_count_cohort_dir,
          showWarnings = FALSE ,
          recursive = T)
        
        # save KM results
        saveRDS(
          object = counts_df, 
          file = file.path(
            out_count_cohort_dir, 
            paste0(
              "counts_", 
              cohorts[j], 
              "_", 
              surv_analysis, 
              ".rds"
            )
          )
        )
      }
      
      # summarize ----
      
      # count number of genes where expression range shows extreme prognosis
      # distinguishing based on distribution shape
      shape_er_l <- lapply(c("normal", "positive_skewed", "bimodal"), function(distr) {
        
        # subset counts keeping only genes with specific distribution shape
        counts_label_df <- counts_df[counts_df$label == distr, ]
        
        # transform each "to_disentangle" hit into 0
        counts_label_df[ counts_label_df == "to_disentangle" ] <- 0
        
        # if more than one row in counts_label_df
        if ( nrow(counts_label_df) > 0 ) {
          
          # loop across shapes
          er_l <- lapply(c("low", "medium", "high"), function(er) {
            
            # expression range bad
            er_bad <- sum(counts_label_df[[paste0("count_", er)]] == -1)
            
            # expression range good
            er_good <- sum(counts_label_df[[paste0("count_", er)]] == 1)
            
            res <- c(
              "bad" = er_bad, 
              "good" = er_good)
            
          })
          
        } else {
          
          er_l <- rep(
              list(setNames(
                object = rep(0, 2), 
                nm = c("bad", "good")))
              , 3)
        }
        
        # name list
        er_l <- setNames(
          object = er_l, 
          nm = c("low", "medium", "high"))
        
        return(er_l)
    
      })
      
      # assign name to list
      shape_er_l <- setNames(
        object = shape_er_l, 
        nm = c("normal", "positive_skewed", "bimodal")
      )
      
      # unlist 
      shape_er_v <- unlist(shape_er_l)
      
      # subsitute "." into "_"
      names(shape_er_v) <- gsub("\\.", "_", names(shape_er_v))
    
      # vector with results
      res <- setNames(
        object = c(
          cohorts[j],
          length(cols_exp_labeled),
          shape_er_v),
        nm = c(
          "cohort", 
          "n_genes_labeled_tested",
          names(shape_er_v))
        )

    } else {
      
      print("No gene passes the correction\n")
      
      # vector with results
      res <- setNames(
        object = c(
          cohorts[j],
          length(cols_exp_labeled), 
          rep(0, 18)), 
        nm = c(
          "cohort",
          "n_genes_labeled_tested",
          "normal_low_bad", 
          "normal_low_good",
          "normal_medium_bad", 
          "normal_medium_good", 
          "normal_high_bad", 
          "normal_high_good", 
          "positive_skewed_low_bad",
          "positive_skewed_low_good", 
          "positive_skewed_medium_bad", 
          "positive_skewed_medium_good",
          "positive_skewed_high_bad", 
          "positive_skewed_high_good", 
          "bimodal_low_bad", 
          "bimodal_low_good", 
          "bimodal_medium_bad", 
          "bimodal_medium_good", 
          "bimodal_high_bad", 
          "bimodal_high_good")
        )
      
    }
    
  } else {
    
    cat(
      "No survival-expression file or labeled genes for", 
      cohorts[j], " is there\n")
    
    res <- NULL
    
  } 
  
  return(res)
  
})


# pool results
count_summary_cohort_df <- as.data.frame(
  do.call(rbind, count_summary_cohort_l))

# transform count column into numeric
count_summary_cohort_df <- count_summary_cohort_df %>%
  mutate_at(
    colnames(count_summary_cohort_df)[-1],
    as.numeric)

# save one summary file for each shape
lapply(c("normal", "positive_skewed", "bimodal"), function(distr) {
  
  # subset count_summary_cohort_df picking only columns of specific shape
  count_summary_cohort_distr_df <- count_summary_cohort_df[
    , c(1, 2, grep(distr, colnames(count_summary_cohort_df)))]
  
  # output file with counts of all cohorts
  output_file <- file.path(
    out_count_dir,
    paste0(
      "count_summary_", 
      distr, 
      "_all-cohorts.rds")
  )
  
  # save
  saveRDS(
    count_summary_cohort_distr_df,
    file = output_file
  )
  
})


# end time
end_time <- Sys.time()