# loop over cohorts
# look at genes where expression range is associated with extreme prognosis
# stratify patients of cohorts in two groups: 
# - "patients_in_the_range_where_geneX_associated_with_extreme"
# - "all_the_other_patients"
# Compare the Ki67 proliferation marker expression in the two stratifications



# import packages ----
library(dplyr)
library(ggplot2)
library(ggpubr)


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

# significance threshold for pairwise KM comparison
kt <- 0.05

# directory with count tables
count_dir <- file.path(
  "../analyses/counts/new/shape-wise",
  paste0(
    as.character(qt),
    "-qvalue-threshold"), 
  paste0(
    as.character(kt), 
    "-km-pairwise-threshold"
  ))

# directory with surv + expression categorical
surv_exp_cat_dir <- file.path(
  "../data/survival-tables/process/surv-time_events-censored_no-imp/survival-exp/categorical_exp/new/shape-wise",
  paste0(
    as.character(qt),
    "-qvalue-threshold"))

# directory with surv + expression continuous
surv_exp_dir <- "../data/survival-tables/process/surv-time_events-censored_no-imp/survival-exp/continuos_exp"


# output variables ----

# directory with info on proliferation status when splitting according
# to genes where expression range is associated with survival
prol_dir <- file.path(
  "../analyses/proliferation/new/shape-wise",
  paste0(
    as.character(qt), 
    "-qvalue-threshold"), 
  paste0(
    as.character(kt), 
    "-km-pairwise-threshold"))

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
  
  # genes_sig pathway
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
  
  # surv_exp_cat pathway 
  surv_exp_cat_path <- file.path(
    surv_exp_cat_dir,
    cohorts[i], 
    surv_analysis, 
    paste0(
      "surv-time_categorical-exp_", 
      surv_analysis, 
      "_",
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
  if ( file.exists(count_path) & 
       file.exists(surv_exp_cat_path) &
       file.exists(surv_exp_cont_path) ) {
    
    # count table
    count_df <- readRDS(count_path)
    
    # surv-exp-categorical df
    surv_exp_cat <- readRDS(surv_exp_cat_path)
    
    # surv-exp df
    surv_exp <- readRDS(surv_exp_cont_path)
    
    
    # wrangle ----
    
    # exclude "to_disentangle" hits in count_low, count_medium and count_high
    # example of why you can find "to_disentangle": 
    # low is significantly different from both medium and high
    # but at that time point its survival is identical to either 
    # medium or high and to disentangle the question we need to plot km.
    # since medium and high are not significant, we do not have
    # to spot which curve has better survival at that time.
    count_sub <- count_df %>% 
      filter(
        count_low != "to_disentangle" , 
        count_medium != "to_disentangle", 
        count_high != "to_disentangle")
    
    ## remove rows where count_low AND count_medium AND count_high == "0"
    # These are the cases where only one comparison is significant
    # for instance low-high comparison is significant, 
    # but low-medium and medium-high comparisons are not significant
    
    # convert last three columns into numeric
    count_sub[, 9:ncol(count_sub)] <- lapply(
      9:ncol(count_sub), function(x) as.numeric(count_sub[[x]])) 
    
    # vector with rows to keep in count_sub
    # that is rows where not all 0 in count_low, count_medium and count_high
    bool_v <- rowSums(abs(count_sub[, 9:ncol(count_sub)])) != 0
    
    # subset count_sub removing rows according to bool_v
    count_sub_sub <- count_sub[bool_v, ]
    
    # subset expression picking only proliferation marker gene 
    # keep only proliferation marker
    surv_exp_sub <- surv_exp[, 
                             colnames(surv_exp) %in% c(
                               "Sample ID", 
                               "TCGA PanCanAtlas Cancer Type Acronym", 
                               "Patient ID", 
                               prol_m)]
    
    # check if samples are in the same order between
    # surv_exp_sub and surv_exp_cat
    all(surv_exp_sub$`Sample ID` == surv_exp_cat$Sample.ID)
    
    # free some space removing surv_exp
    rm(surv_exp)
    
    
    # loop ----
    
    # loop over genes in count
    agreement_info_l <- lapply(seq(dim(count_sub_sub)[1]), function(j) {
      
      cat("Process", count_sub_sub$gene[j], "\n")
      
      # check in which range the gene is associated with extreme prognosis
      # start with "low" expression range
      # by running the following snippet whenever both "low" and "high"
      # or "low and "medium" have an effect on survival, I will always compare
      # Ki67 expression distributions stratifying patients based on "low" group. 
      # I should get a mirrored result if instead I stratify based on "medium" or "high" group, 
      # but I cannot be sure about that
      if ( as.numeric(count_sub_sub[j, "count_low"]) != 0 ){ # if there is a prognostic effect in low
        
        if ( as.numeric(count_sub_sub[j, "count_low"]) == 1 ) { # if the prognostic effect is good
          
          info <- c(
            "cohort" = cohorts[i],
            "gene" = count_sub_sub$gene[j], 
            "label" = count_sub_sub$label[j], 
            "range" = "low", 
            "prognostic_effect" = "good")
          
        } else { # the prognostic effect is bad
          
          info <- c(
            "cohort" = cohorts[i],
            "gene" = count_sub_sub$gene[j], 
            "label" = count_sub_sub$label[j], 
            "range" = "low",  
            "prognostic_effect" = "bad")
        }
        
      } else if ( as.numeric(count_sub_sub[j, "count_medium"]) != 0 ) { # is there is a prognostic effect in medium
        
        if ( as.numeric(count_sub_sub[j, "count_medium"]) == 1 ) {
          
          info <- c(
            "cohort" = cohorts[i],
            "gene" = count_sub_sub$gene[j], 
            "label" = count_sub_sub$label[j], 
            "range" = "medium", 
            "prognostic_effect" = "good")
          
        } else {
          
          info <- c(
            "cohort" = cohorts[i],
            "gene" = count_sub_sub$gene[j], 
            "label" = count_sub_sub$label[j], 
            "range" = "medium", 
            "prognostic_effect" = "bad")
        }
        
      } else if ( as.numeric(count_sub_sub[j, "count_high"]) != 0 ) {
        
        if ( as.numeric(count_sub_sub[j, "count_high"]) == 1 ) {
          
          info <- c(
            "cohort" = cohorts[i],
            "gene" = count_sub_sub$gene[j], 
            "label" = count_sub_sub$label[j], 
            "range" = "high", 
            "prognostic_effect" = "good")
          
        } else {
          
          info <- c(
            "cohort" = cohorts[i],
            "gene" = count_sub_sub$gene[j], 
            "label" = count_sub_sub$label[j], 
            "range" = "high", 
            "prognostic_effect" = "bad")
        }
      }
      
      # patients stratification based on gene j
      strata <- as.character(
        surv_exp_cat[[grep(count_sub_sub$gene[j], colnames(surv_exp_cat))]])
      
      # transform stratification from 3 to 2 groups 
      if ( info[names(info) == "range"] == "low" ) {
        
        # transform what is not low in "non-low"
        strata_transformed <- factor(
          ifelse(
            strata == "low", "low", "non-low"),
          levels = c("low", "non-low"))
        
      } else if ( info[names(info) == "range"] == "medium") {
        
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
      
      # add column with stratification
      surv_exp_sub$strata <- strata_transformed
      
      ## compare means of Ki67 expression distributions
      res <- compare_means(
        formula = as.formula(paste0(prol_m, " ~ strata")),
        data = surv_exp_sub)
      
      # if the comparison is significant
      # compute means 
      if (res$p.adj < 0.05) {
        
        if ( info[names(info) == "range"] == "low" ) {
          
          # mean group1
          mean_group1 <- mean(surv_exp_sub[[prol_m]][surv_exp_sub$strata == "low"])
          
          # mean group2
          mean_group2 <- mean(surv_exp_sub[[prol_m]][surv_exp_sub$strata == "non-low"])
          
        } else if ( info[names(info) == "range"] == "medium" ) {
          
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
          if (info[names(info) == "prognostic_effect"] == "good") {
            
            agreement_surv_prol <- "agree"
            
          } else {
            agreement_surv_prol <- "disagree"
          }
        } else { # mean(group1) > mean(group2)
          
          if (info[names(info) == "prognostic_effect"] == "bad") {
            agreement_surv_prol <- "agree"
            
          } else {
            agreement_surv_prol <- "disagree"
          }
        }
        
        # add agreement_surv_prol to info
        info <- c(
          info, 
          "agreement_surv_prol" = agreement_surv_prol
        )
        
      } else {
        
        cat(
          names(prol_m),
          "proliferation does not significantly change between the two groups\n")
        
        # add agreement_surv_prol to info
        info <- c(
          info, 
          "agreement_surv_prol" = paste0(
            names(prol_m), 
            "_proliferation_not_significant")
        )
      }
    })
    
    # rbind
    agreement_info_df <- as.data.frame(do.call(rbind, agreement_info_l))
    
    
  } else {
    
    cat("Needed files are not there\n")
    agreement_info_df <- NULL
  }
  
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
      "-proliferation_genes_expression-range-associated-with-extreme-prognosis.rds")
  ))