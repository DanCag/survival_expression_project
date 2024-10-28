# test association between bimodal distribution shape 
# and extreme prognosis in medium group
# contingency table with bimodal 
#
#            extreme prognosis in medium    extreme prognosis not in medium 
# bimodal
# non-bimodal 


# number of samples info is in "../data/survival-tables/info/info_extracted-samples_DSS"

# packages ----
library(dplyr)


# parameters ----

# survival analysis
surv_analysis <- "DSS"

# significance threshold for multiple testing correction 
qvalue_threshold <- 0.2

# significance threshold for pairwise KM comparison
km_pairwise_threshold <- 0.05

# method for multiple test correction
method <- "BH"

# split by
# it could be
# - "tertiles"
# - "shape-wise"
split_by <- "shape-wise"


# input ----

# count directory
count_dir <- file.path(
  "../analyses/counts/new",
  split_by, 
  paste0(
    as.character(qvalue_threshold),
    "-qvalue-threshold"), 
  paste0(
    as.character(km_pairwise_threshold), 
    "-km-pairwise-threshold"))


# output ----

# output directory
out_dir <- file.path(
  "../output/analyses/kaplan-meier/association/counts/new", 
  split_by)

# make output directory
dir.create(
  out_dir, 
  showWarnings = F, 
  recursive = T)


# data ----

# cancer cohorts
cohorts <- readRDS("../output/cohorts.rds")


# loop ----

# loop over cohorts and run Fisher test
fisher_cohort_l <- lapply(seq(length(cohorts)), function(i) {
  
  cat("Process", cohorts[i], "cohort\n", 
      "--------------------", 
      "\n\n")
  
  # count file
  count_path <- file.path(
    count_dir, 
    cohorts[i], 
    surv_analysis,
    paste0(
      "counts_", 
      cohorts[i],
      "_", 
      surv_analysis, 
      ".rds")
  )
  
  # check if count file exists
  if ( file.exists(count_path) ) {
    
    # import count file
    count_df <- readRDS(count_path)
    
    # wrangling ----
    
    # remove "to_disentangle" from count_df
    count_sub_df <- count_df %>% 
      filter(
        count_low != "to_disentangle", 
        count_medium != "to_disentangle",
        count_high != "to_disentangle")
    
    # transform count columns into numeric
    count_sub_df <- count_sub_df %>% 
      mutate_at(
        c("count_low", "count_medium", "count_high"),
        as.numeric)
    
    # remove genes with 0 in count_low, count_medium and count_high
    # these are genes for which only one comparison is significant, like 
    # for instance high and low, but medium-low and medium-high are not significant, 
    # therefore we cannot really spot an extreme prognostic behavior
    to_keep <- rowSums(abs(count_sub_df[, 9:11])) != 0
    count_sub_df <- count_sub_df[to_keep, ]
    
    
    # Fisher's ----
    
    ## shape medium associated
    
    # genes with shape x associated with survival in medium
    shape_medium_ass_genes <- count_sub_df$gene[
      count_sub_df$label == "bimodal" & 
        count_sub_df$count_medium %in% c(-1, 1)]
    
    ## shape medium not associated
    
    # number of genes with specific shape and medium not associated with survival
    shape_medium_not_ass_genes <- count_sub_df$gene[
      count_sub_df$label == "bimodal" & 
        count_sub_df$count_medium == 0]
    
    ## not shape associated
    
    # genes with shape different from x associated with survival
    not_shape_medium_ass_genes <- count_sub_df$gene[
      count_sub_df$label != "bimodal" & 
        count_sub_df$count_medium %in% c(-1, 1)]
    
    ## not shape not medium ass
    
    # genes with shape different from x medium not associated with survival
    not_shape_medium_not_ass_genes <- count_sub_df$gene[
      count_sub_df$label != "bimodal" & 
        count_sub_df$count_medium == 0]
    
    # contigency table
    m <- matrix(
      c(length(shape_medium_ass_genes), 
        length(shape_medium_not_ass_genes), 
        length(not_shape_medium_ass_genes), 
        length(not_shape_medium_not_ass_genes)), 
      nrow = 2, 
      ncol = 2, 
      byrow = T, 
      dimnames = list(
        c("bimodal", "not_bimodal"),
        c("medium_associated", "medium_not_associated")))
    
    fisher <- fisher.test(m)
    pval <- fisher$p.value
    estimate <- fisher$estimate
    
    # chisquare
    # chisq_res <- chisq.test(m)
    
    # combine pvalue and estimate
    res <- c(
      "cohort" = cohorts[i],
      "pvalue" = pval,
      estimate)
    
  } else {
    
    cat("No files for ", cohorts[i], "\n")
    
    res <- c(
      "cohort" = cohorts[i],
      rep(NA, 2))
    
  }
  
  return(res)
})

# rbind elements of list
fisher_cohort_df <- as.data.frame(do.call(rbind, fisher_cohort_l))

# transform columns into numeric
fisher_cohort_df <- fisher_cohort_df %>% 
  mutate_at(c(2:3), as.numeric)

# change columns' names
colnames(fisher_cohort_df) <- sub(" ", "_", colnames(fisher_cohort_df))

# save table with associations
write.table(
  fisher_cohort_df, 
  file = file.path(
    out_dir, 
    paste0(
      "association_shape-medium-extreme_", 
      split_by, 
      ".tsv"
    )
  ), 
  sep = "\t", 
  row.names = F, 
  col.names = T
)
