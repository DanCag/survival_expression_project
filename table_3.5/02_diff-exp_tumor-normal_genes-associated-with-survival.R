# compare expression of genes associated with survival
# between tumor and normal samples in those cohorts
# where there are at least 20 normal samples


# packages ----
library(dplyr)
library(ggplot2)
source("./functions/translate_prognosis.R")
source("./functions/classification_f.R")


# parameters ----

# significance threshold for multiple testing correction
qt <- 0.2

# method for multiple test correction
correction_method <- "BH"

# significance threshold for log-rank test pairwise comparisons
kmpt <- 0.05

# survival analysis
surv_analysis <- "DSS"

# paired or unpaired?
# choose how to run wilcoxon test
wilcoxon_strategy <- "paired"


# input ----

# file with info about the number of normal samples in each cohort
normal_info <- read.delim(
  "../data/expression_normal_samples/tcga_normal/normal_samples_info.tsv")

# directory with count 
count_dir <- file.path(
  "../analyses/counts/new/shape-wise",
  paste0(
    as.character(qt),
    "-qvalue-threshold"),
  paste0(
    as.character(kmpt),
    "-km-pairwise-threshold"))

# tumor expression directory
surv_exp_dir <- "../data/survival-tables/process/surv-time_events-censored_no-imp/survival-exp/continuos_exp"

# normal expression directory
normal_exp_dir <- "../data/expression_normal_samples/tcga_normal"


# output ----

# output directory
out_dir <- file.path(
  "../analyses/comparison_tumor-normal/quantitative/significant", 
  wilcoxon_strategy)


# wrangle ----

# pick only those cohorts where the number of normal samples > 20
cohorts <- normal_info$cohort[normal_info$n_normal_samples > 0 & normal_info$n_normal_samples <= 20]


# loop ----

# loop over cohorts and compare tumor and normal expression for genes 
# with an expression range associated with extreme prognosis
gene_comparison_l <- lapply(seq(cohorts), function(j) {
  
  # loop over cohorts
  cat("Process", cohorts[j], "cohort\n")
  
  # count filepath
  count_path <- file.path(
    count_dir,
    cohorts[j], 
    surv_analysis, 
    paste0(
      "counts_",
      cohorts[j],
      "_", 
      surv_analysis, 
      ".rds"))
  
  # tumor surv-exp filepath
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
  
  # normal exp
  normal_exp_path <- file.path(
    normal_exp_dir, 
    cohorts[j],
    paste0(
      "tcga_normal_", 
      cohorts[j], 
      ".rds")
  )
  
  # later add pvalue upon randomization table
  
  # check if surv-exp file exists for that specific analysis
  if ( file.exists(count_path) &
       file.exists(surv_exp_path) &
       file.exists(normal_exp_path) ) {
    
    # import datasets ----
    
    # count
    count <- readRDS(count_path)
    
    # tumor survival + expression
    surv_exp <- readRDS(surv_exp_path)
    
    # normal expression
    normal_exp <- readRDS(normal_exp_path)
    
    
    # wrangle ----
    
    # exclude "to_disentangle" hits 
    count_sub <- count %>% 
      filter(
        count_low != "to_disentangle", 
        count_medium != "to_disentangle", 
        count_high != "to_disentangle")
    
    # genes to remove because they don't have any effect on survival
    to_keep <- rowSums(abs(count_sub[, c(9, 10, 11)])) != 0
    
    # subset count
    count_sub_sub <- count_sub[to_keep, ]
    
    # count_sub_sub must not be empty
    if ( nrow(count_sub_sub) > 0 ) {
      
      # add column with patient to normal_exp
      normal_exp <- normal_exp %>% 
        mutate(patient_id = substr(sample_submitter_id, 1, 12), 
               .after = sample_submitter_id)
      
      # if paired comparison
      if (wilcoxon_strategy == "paired") {
        
        # patients with both tumor and normal samples
        match_patients <- intersect(surv_exp$`Patient ID`, normal_exp$patient_id)
        
        # order
        match_patients <- match_patients[order(match_patients)]
        
        if (length(match_patients) != 0 ){
          
          # subset tumor
          tumor_df <- surv_exp[
            surv_exp$`Patient ID` %in% match_patients, c(1, 8:ncol(surv_exp))]
          
          # subset normal
          normal_df <- normal_exp[
            normal_exp$patient_id %in% match_patients, ]
          
        } else {
          
          tumor_df <- NULL
          normal_df <- NULL
        }
        
      } else {
        
        # tumor
        tumor_df <- surv_exp[, c(1, 8:ncol(surv_exp))]
        
        #  normal
        normal_df <- normal_exp
      }
      
      
      if (!is.null(tumor_df) & !is.null(normal_df) ) {
        
        # loop over genes
        gene_comparison_l <- lapply(count_sub_sub$gene, function(g) {
          
          # cat("gene:", g, "\n")
          
          # translate prognosis into a string
          prognosis <- translate_prognosis(count_sub_sub, gene = g)
          
          # dataframe with tumor and normal expression
          df <- data.frame(
            "exp" = c(tumor_df[[g]], normal_df[[g]]), 
            "sample_type" = factor(
              c(
                rep("tumor", length(tumor_df[[g]])), 
                rep("normal", length(normal_df[[g]])))
              , 
              c(levels = "tumor", "normal")))
          
          # check if df is all 0
          if ( all(df$exp == 0) ) {
            
            # cat("All tumor and normal samples have 0 expression\n")
            
            res <- c(
              "gene" = g, 
              "tumor_normal" = NA,
              "pvalue" = NA)
            
          } else {
            
            # compare expression ----
            
            if (wilcoxon_strategy == "paired") {
              
              # wilcoxon test
              # wilcox_res <- compare_means(
              #   formula = exp ~ sample_type, 
              #   data = df, 
              #   method = "wilcox.test")
              wilcox_res <- wilcox.test(
                x = df$exp[df$sample_type == "tumor"],
                y = df$exp[df$sample_type == "normal"], 
                paired = T)
              
            } else {
              
              wilcox_res <- wilcox.test(
                x = df$exp[df$sample_type == "tumor"],
                y = df$exp[df$sample_type == "normal"], 
                paired = F)
            }
            
            
            # if the comparison is significant
            # compute means 
            if (wilcox_res$p.value < 0.05) {
              
              # mean tumor
              mean_tumor <- mean(df$exp[df$sample_type == "tumor"])
              
              # mean normal
              mean_normal <- mean(df$exp[df$sample_type == "normal"])
              
              if (mean_tumor > mean_normal) {
                
                tumor_normal <- "higher"
                
              } else if (mean_tumor < mean_normal) {
                
                tumor_normal <- "lower"
                
              } else { # should not occur
                
                tumor_normal <- "no_difference"
                
              }
              
            } else {
              
              tumor_normal <- "no_difference"
            }
            
            # result
            res <- c(
              "gene" = g, 
              "prognosis_all_tumor_samples" = prognosis, 
              "tumor_normal" = tumor_normal,
              "pvalue" = wilcox_res$p.value)
            
          }
        })
        
        # rbind elements
        gene_comparison_df <- bind_rows(gene_comparison_l)
        
        # pvalue correction 
        gene_comparison_df <- gene_comparison_df %>% 
          mutate(
            p_adj = p.adjust(
              pvalue, 
              method = correction_method), 
            .after = pvalue)
        
        
        # save ----
        
        # directory where to store iteration result
        out_cohort_dir <- file.path(
          out_dir, 
          cohorts[j])
        
        # create directory(ies) recursively if not there
        dir.create(
          out_cohort_dir,
          showWarnings = FALSE ,
          recursive = T)
        
        # save res 
        saveRDS(
          object = gene_comparison_df, 
          file = file.path(
            out_cohort_dir, 
            paste0(
              "comparison_significant_tumor-normal_",
              wilcoxon_strategy, 
              "_",
              cohorts[j],
              ".rds")
          )
        )
        
      } else {
        
        cat("No patients in tumor_df and normal_df\n")
      }
      
      
    } else {
      
      cat("No significant genes for", cohorts[j], "\n")
    }
    
  } else {
    
    cat(
      "Needed files for cohort", 
      cohorts[j], 
      "are not there\n")
    
    # vector with results
    gene_comparison_df <- NULL
    
  }
  
})
