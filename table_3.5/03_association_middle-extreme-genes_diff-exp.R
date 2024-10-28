# test association between extreme prognosis
# in medium group and differential expression
# example of contingency table 
#
#                      diff. expression YES  diff. expression NO
# middle-extreme YES
# non-significant



# packages ----
library(dplyr)


# parameters ----

# significance threshold for multiple testing correction
qt <- 0.2

# paired or unpaired?
# choose how to run wilcoxon test: "paired" or "unpaired"
wilcoxon_strategy <- "unpaired"


# input ----

# directory with differential expression info
diff_exp_dir <- "../output/analyses/kaplan-meier/comparison_tumor-normal/quantitative"


# output ----

# output directory
out_dir <- "../output/analyses/kaplan-meier/association/tumor_normal"

# make output directory
dir.create(
  out_dir, 
  showWarnings = F, 
  recursive = T)


# data ----

# file with info about the number of normal samples in each cohort
normal_info <- read.delim(
  "../output/expression/tcga_normal/normal_samples_info.tsv")


# wrangle ----

# pick only those cohorts where the number of normal samples > 20
# cohorts <- normal_info$cohort[normal_info$n_normal_samples > 20]
cohorts <- normal_info$cohort


# loop ----

# loop over cohorts and run Fisher test
fisher_cohort_l <- lapply(seq(length(cohorts)), function(i) {
  
  cat("Process", cohorts[i], "cohort\n", 
      "--------------------", 
      "\n\n")
  
  # significant genes
  sig_path <- file.path(
    diff_exp_dir,
    "significant",
    wilcoxon_strategy, 
    cohorts[i], 
    paste0(
      "comparison_significant_tumor-normal_", 
      wilcoxon_strategy, 
      "_", 
      cohorts[i],
      ".rds")
  )
  
  # non-significant genes
  non_sig_path <- file.path(
    diff_exp_dir,
    "non-significant",
    wilcoxon_strategy, 
    cohorts[i], 
    paste0(
      "comparison_non-significant_tumor-normal_",
      wilcoxon_strategy, 
      "_", 
      cohorts[i],
      ".rds")
  )
  
  # check if files exist
  if ( file.exists(sig_path) &
       file.exists(non_sig_path) ) {
    
    # import significant file
    sig <- readRDS(sig_path)
    
    # import non significant file
    non_sig <- readRDS(non_sig_path)
    
    
    # wrangling ----
    
    # consider only genes 
    # where medium expression range is associated with extreme prog
    sig_sub <- sig %>% 
      filter(grepl("medium", prognosis_all_tumor_samples))
    
    
    # Fisher's ----
    
    # medium-extreme genes differential expressed
    medium_ass_diff_genes <- sig_sub %>% 
      filter(p_adj < qt, tumor_normal %in% c("lower", "higher")) %>% 
      pull(gene)
    
    # medium-extreme genes not differential expressed
    medium_ass_not_diff_genes <- sig_sub %>% 
      filter(tumor_normal == "no_difference") %>% 
      pull(gene)
    
    # non-significant genes differential expressed
    non_sig_diff_genes <- non_sig %>% 
      filter(p_adj < qt, tumor_normal %in% c("lower", "higher")) %>% 
      pull(gene)
    
    # non-significant genes not differential expressed
    non_sig_not_diff_genes <- non_sig %>% 
      filter(tumor_normal == "no_difference") %>% 
      pull(gene)
    
    # contingency table
    m <- matrix(
      c(length(medium_ass_diff_genes), 
        length(medium_ass_not_diff_genes), 
        length(non_sig_diff_genes), 
        length(non_sig_not_diff_genes)), 
      nrow = 2, 
      ncol = 2, 
      byrow = T, 
      dimnames = list(
        c("medium-extreme", "non-significant"),
        c("diff-exp", "not-diff-exp")))
    
    fisher <- fisher.test(m)
    pval <- fisher$p.value
    estimate <- fisher$estimate
    
    # chisquare
    # chisq_res <- chisq.test(m)
    
    # combine pvalue and estimate
    res <- c(
      pval,
      estimate)
    
    # # in case of chisquare
    # res <- c(
    #   "pvalue" = chisq_res$p.value)
    
    # info with cohort and fisher's test results
    info <- c("cohort" = cohorts[i], res)
    
  } else {
    
    cat("No files for ", cohorts[i],"\n")
    
    info <- c(
      "cohort" = cohorts[i],
      rep(NA, 2))
    
    # in case of chisquare
    # info <- c(
    #   "cohort" = cohorts[i], 
    #   rep(NA, 3))
  }
  
  return(info)
})

# rbind elements of list
fisher_cohort_df <- as.data.frame(do.call(rbind, fisher_cohort_l))
colnames(fisher_cohort_df)[2:3] <- c("pvalue", "odds_ratio")

# transform columns into numeric
fisher_cohort_df <- fisher_cohort_df %>% 
  mutate_at(c(2:3), as.numeric)

# save table with associations
write.table(
  fisher_cohort_df, 
  file = file.path(
    out_dir, 
    paste0(
      "association_medium-extreme_diff-exp_",
      wilcoxon_strategy,
      ".tsv")), 
  sep = "\t", 
  row.names = F, 
  col.names = T)
