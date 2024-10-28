# loop over cohorts
# look at genes we could split into "low", "medium", "high" expression range, 
# but where no expression range is associated with extreme prognosis
# (KM analysis failed)
# Stratify patients of cohorts in two groups: 
# - "patients_in_the_range_of_interest" (for instance medium)
# - "all_the_other_patients"
# Compare the Ki67 proliferation marker expression in the two stratifications
# Does the Ki67 expression significantly change in the two groups?


# import packages ----
library(parallel)
library(dplyr)
library(ggplot2)
library(ggpubr)
source("./functions/split_function.R")


# input variables ----

# survival analysis
surv_analysis <- "DSS"

# cancer cohorts
cohorts <- readRDS("../data/cohorts.rds")

# proliferation markers
# MCM proteins family: MCM2-MCM7
prol_m <- c(
  "Ki67" = "ENSG00000148773")

# significance threshold for multiple testing correction 
qt <- 0.2

# method for multiple test correction
method <- "BH"

# range of interest
range_of_interest <- "high"

# split by
split_by <- "quartiles"

# number of cores for parallelization
n_cores <- 20

# directory with non significant genes
not_sig_dir <- file.path(
  "../analyses/genes_list/new/shape-wise/non-significant",
  paste0(
    as.character(qt),
    "-qvalue-threshold"))

# classification info for all genes in cohort
classification_dir <- "../analyses/classification/TCGA"

# directory with surv + expression continuous
surv_exp_dir <- "../data/survival-tables/process/surv-time_events-censored_no-imp/survival-exp/continuos_exp"


# output variables ----

# directory with info on proliferation status when splitting according
# to genes where expression range is associated with survival
prol_dir <- "../analyses/proliferation/new/shape-wise/non-significant"

# create directory(ies) recursively if not there
dir.create(
  prol_dir,
  showWarnings = FALSE ,
  recursive = T)


# loop ----

# loop over cohorts
agreement_info_cohort_l <- lapply(seq(cohorts), function(i) {
  
  cat("Process", cohorts[i], "cohort\n")
  
  # import datasets ----
  
  # genes_non_sig pathway
  not_sig_path <- file.path(
    not_sig_dir,
    cohorts[i], 
    surv_analysis, 
    paste0(
      "genes_non-significant-after-mtc_", 
      cohorts[i],
      "_",
      surv_analysis, 
      ".txt"))
  
  # classification path
  classification_path <- file.path(
    classification_dir, 
    cohorts[i], 
    paste0(
      "genes_classification_", 
      cohorts[i], 
      ".rds"))
  
  # surv_exp_continuous pathway (need to get proliferation marker expression)
  surv_exp_cont_path <- file.path(
    surv_exp_dir,
    cohorts[i], 
    surv_analysis, 
    paste0(
      "surv-time_exp_events-censored_no-imp_", 
      surv_analysis, 
      "_",
      cohorts[i], 
      ".rds"))
  
  # check if all the three files exist
  if ( file.exists(not_sig_path) & 
       file.exists(classification_path) &
       file.exists(surv_exp_cont_path) ) {
    
    # non signficant genes
    not_sig_genes <- read.csv(not_sig_path)[[1]]
    
    cat("Non significant genes to parse:", length(not_sig_genes), "\n")
    
    # surv-exp df
    surv_exp <- readRDS(surv_exp_cont_path)
    
    # classification info
    classification_info <- readRDS(classification_path)
    
    # subset expression picking only proliferation marker gene
    # keep only proliferation marker
    surv_exp_sub <- surv_exp[,
                             colnames(surv_exp) %in% c(
                               "Sample ID",
                               "TCGA PanCanAtlas Cancer Type Acronym",
                               "Patient ID",
                               prol_m)]
    
    # not_sig_genes <- not_sig_genes[1:10]
    
    
    # loop ----
    
    # loop over genes in count
    agreement_info_l <- mclapply(seq(length(not_sig_genes)), function(j) {
      
      # cat("Process", not_sig_genes[j], "\n")
      
      # shape of gene
      shape <- classification_info$label[
        classification_info$gene == not_sig_genes[j]]
      
      if ( shape %in% c("normal", "positive_skewed") ) {
        
        # Know from 01_split.R that this gene can be split in four groups 
        
        # bin expression in four groups
        strata <- split_f(
          gene = not_sig_genes[j], 
          split_strategy = split_by, 
          surv_exp_df = surv_exp)
        
      } else { 
        
        # we deal with bimodal genes
        # know already from 01_split.R that this gene can be split in 3 groups
        # and that in each group there are minimum 3 samples
        
        # bin expression in three groups
        strata <- split_bimodal3_f(
          gene = not_sig_genes[j],
          surv_exp_df = surv_exp)
        
      }
      
      # collect info
      info <- c(
        "cohort" = cohorts[i],
        "gene" = not_sig_genes[j], 
        "label" = shape, 
        "range" = "none", 
        "prognostic_effect" = "none") 
      
      # transform stratification from 3 to 2 groups 
      if ( range_of_interest == "low" ) {
        
        # transform what is not low in "non-low"
        strata_transformed <- factor(
          ifelse(
            strata == "low", "low", "non-low"),
          levels = c("low", "non-low"))
        
      } else if ( range_of_interest == "medium" ) {
        
        # transform what is not medium in "non-medium"
        strata_transformed <- factor(
          ifelse(
            strata == "medium", "medium", "non-medium"),
          levels = c("medium", "non-medium"))
        
      } else {
        
        # transform what is not high in "non-high"
        strata_transformed <- factor(
          ifelse(
            strata == "high", "high", "non-high"),
          levels = c("high", "non-high"))
        
      }
      
      # to surv_exp_sub add column with stratification
      surv_exp_sub$strata <- strata_transformed
      
      ## compare means of Ki67 expression distributions
      res <- compare_means(
        formula = as.formula(paste0(prol_m, " ~ strata")),
        data = surv_exp_sub)
      
      # if the comparison is significant
      # compute means 
      if (res$p.adj < 0.05) {
        
        if ( range_of_interest == "low" ) {
          
          # mean group1
          mean_group1 <- mean(surv_exp_sub[[prol_m]][surv_exp_sub$strata == "low"])
          
          # mean group2
          mean_group2 <- mean(surv_exp_sub[[prol_m]][surv_exp_sub$strata == "non-low"])
          
        } else if ( range_of_interest == "medium" ) {
          
          # mean group1
          mean_group1 <- mean(surv_exp_sub[[prol_m]][surv_exp_sub$strata == "medium"])
          
          # mean group2
          mean_group2 <- mean(surv_exp_sub[[prol_m]][surv_exp_sub$strata == "non-medium"])
          
        } else {
          
          # mean group1
          mean_group1 <- mean(surv_exp_sub[[prol_m]][surv_exp_sub$strata == "high"])
          
          # mean group2
          mean_group2 <- mean(surv_exp_sub[[prol_m]][surv_exp_sub$strata == "non-high"])
          
        }
        
        # compare means
        if ( mean_group1 < mean_group2 ) {
          
          # prognostic effect
          effect <- paste0(range_of_interest, "_less_than_non_", range_of_interest)
          
        } else { # mean(group1) > mean(group2)
          
          # prognostic effect
          effect <- paste0(range_of_interest, "_more_than_non_", range_of_interest)
        }
        
        # add agreement_surv_prol to info
        info <- c(
          info, 
          "proliferation_change" = "yes",
          "proliferation_effect" = effect)
        
      } else {
        
        # cat(
        #   names(prol_m),
        #   "proliferation does not significantly change between the two groups\n")
        
        # add agreement_surv_prol to info
        info <- c(
          info, 
          "proliferation_change" = "no",
          "proliferation_effect" = "no")
      }
    }, mc.cores = n_cores)
    
    # free some space removing surv_exp
    rm(surv_exp)
    
    # rbind
    agreement_info_df <- as.data.frame(do.call(rbind, agreement_info_l))
    
    
  } else {
    
    cat("Needed files are not there\n")
    agreement_info_df <- NULL
  }
  
  cat("\n")
  
  # return agreement ino_df
  return(agreement_info_df)
  
})

# remove null elements from list
agreement_info_cohort_no_null_l <- agreement_info_cohort_l[
  -which(sapply(agreement_info_cohort_l, is.null))]

# pool dataframes into one
agreement_info_cohort_df <- do.call(rbind, agreement_info_cohort_no_null_l)

# save it
saveRDS(
  agreement_info_cohort_df,
  file.path(
    prol_dir,
    paste0(
      "agreement_survival-",
      names(prol_m),
      "-proliferation_referred-to-", 
      range_of_interest, 
      "_non-significant-genes.rds")
  ))
