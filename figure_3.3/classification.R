# loop over cancer cohorts and classify gene expression distributions into: 
# - bimodal
# - normal
# - skewed (positively or negatively)

# according to the shape, split samples in a specific way
# get a table with the number of genes per each classification group
# use it to make figure 3.3 of thesis

# run this script on workstation

# import packages ----
library(parallel)
source("./functions/classification_f.R")
library(ggplot2)
library(gridExtra)
library(dplyr)


# input variables ----

# cancer cohorts
cohorts <- readRDS("../data/cohorts.rds")

# survival + expression (continuous) table
surv_exp_cont <- "../data/survival-tables/process/surv-time_events-censored_no-imp/survival-exp/continuos_exp"

# survival analysis
surv_analysis <- "DSS"

# number of cores
n_cores <- 20


# output variables ----

# pathway where we will store info about genes tested
info_genes_path <- "../analyses/classification/TCGA/info_shape_genes.tsv"


# loop ----

# loop over cohort 
genes_df_l <- lapply(seq(cohorts), function(i) {
  
  cat("\n\n", "Process", cohorts[i], "cohort","\n\n", 
      "------------------", 
      "\n\n")
  
  # surv-exp filepath
  surv_exp_path <- file.path(
    surv_exp_cont,
    cohorts[i], 
    surv_analysis, 
    paste0(
      "surv-time_exp_events-censored_no-imp_", 
      surv_analysis, 
      "_",
      cohorts[i], 
      ".rds")
  )
  
  # check if surv-exp file exists for that specific analysis
  # if expression matrix exists
  if ( file.exists(surv_exp_path) ) {
    
    # import datasets ----
    
    # survival + expression
    surv_exp <- readRDS(surv_exp_path)
    
    
    # Data wrangling ----
    
    # subset cols with gene expression
    cols_exp <- colnames(surv_exp)[
      grep("ENS", colnames(surv_exp))]
    
    # remove genes with all 0s 
    cols_exp_selected <- cols_exp[
      apply(surv_exp[, cols_exp], 
            2, 
            function(x) sum(x) != 0)]
    
    
    # loop ----
    
    # loop over genes and label them
    distr_classification_l <- mclapply(
      cols_exp_selected,
      function(g) {
        
        # cat("Assess", g, "gene distribution\n")
        cl <- classification_f(gene = g, exp = surv_exp)
        
      }, mc.cores = n_cores)
    
    # dataframe with gene, label and pvalue information
    label_df <- data.frame(
      "gene" = sapply(distr_classification_l, "[[", 1), 
      "label" = sapply(distr_classification_l, "[[", 2))
    

    ## genes info ----
    
    # discard unlabeled genes
    label_df_no_unlabeled_df <- label_df[label_df$label != "unlabeled", ]
    
    # genes info dataframe
    genes_info_df <- data.frame(
      "n_starting_genes" = length(cols_exp),
      "n_tested_genes" = length(cols_exp_selected),
      "n_labeled_genes" = dim(label_df_no_unlabeled_df)[1],
      rbind(table(label_df$label))
    )
    
    # rename last three columns
    colnames(genes_info_df)[4:ncol(genes_info_df)] <- paste0(
      "n_",
      colnames(genes_info_df)[4:ncol(genes_info_df)],
      "_genes")

  } else {
    
    cat("No survival-expression file for that analysis is there\n")
    genes_info_df <- NULL
  }
  
  return(genes_info_df)
})


# info on genes ----

# assign cohorts names to genes_df_l
genes_df_l <- setNames(
  object = genes_df_l,
  nm = cohorts)

# rbind genes_df_l elements in one dataframe
genes_info_df <- bind_rows(genes_df_l)

# index of cohort to discard (LAML) because NULL element of list 
idx_to_discard <- which(
  cohorts == names(genes_df_l)[unlist(lapply(genes_df_l, is.null))])

# add cohorts info
genes_info_df <- data.frame(
  "cohort" = cohorts[-idx_to_discard],
  genes_info_df)

# replace NAs with 0
genes_info_df[is.na(genes_info_df)] <- 0

# save it
write.table(
  genes_info_df,
  file = info_genes_path,
  sep = "\t",
  row.names = F,
  col.names = T)


# summarise ----

# compute fraction of classification classes
cl_freq <- genes_info_df %>% 
  mutate(
    freq_bimodal = n_bimodal_genes/n_tested_genes, 
    freq_normal = n_normal_genes/n_tested_genes, 
    freq_positive_skewed = n_positive_skewed_genes/n_tested_genes, 
    freq_negative_skewed = n_negative_skewed_genes/n_tested_genes, 
    freq_unlabeled_skewed = n_unlabeled_genes/n_tested_genes
  ) %>% 
  select(
    cohort, 
    freq_bimodal,
    freq_normal, 
    freq_positive_skewed, 
    freq_negative_skewed, 
    freq_unlabeled_skewed)

# melt 
cl_freq_gathered <- cl_freq %>% 
  gather("shape", "frequency", -cohort)


# barplot ----

# frequency barplot
ggplot(data = cl_freq_gathered) +
  geom_bar(
    mapping = aes(x = forcats::fct_rev(cohort), y = frequency, fill = shape), 
    stat = "identity",
    position = "stack") +
  scale_fill_manual(values = c("#ffbc42", "#d81159", "#8f2d56", "#218380", "#73d2de")) + 
  labs(y = "frequency", x = "cohort") +
  theme_bw() +
  coord_flip() +
  geom_text(
    data = genes_info_df,
    aes(x = cohort, y = 0, label = n_tested_genes))

