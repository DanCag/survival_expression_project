# association survival-proliferation


# packages ----
library(dplyr)

# parameters ----

# significance threshold for multiple testing correction 
qt <- 0.2

# significance threshold for KM analysis 
kt <- 0.05

# range of interest
range_of_interest <- "high"

# output directory
out_dir <- "../analyses/association/proliferation/new"


# data ----

# agreement significant genes
agreement_sig <- readRDS(
  file.path(
    "../analyses/proliferation/new/shape-wise", 
    paste0(
      as.character(qt),
      "-qvalue-threshold"), 
    paste0(
      as.character(kt),
      "-km-pairwise-threshold"), 
    "agreement_survival-Ki67-proliferation_genes_expression-range-associated-with-extreme-prognosis.rds"
  )
)

# agreement non-significant
agreement_non_sig <- readRDS(
  file.path(
    "../analyses/proliferation/new/shape-wise/non-significant", 
    paste0(
      "agreement_survival-Ki67-proliferation_referred-to-",
      range_of_interest,
      "_non-significant-genes.rds")))

# cancer cohorts
cohorts <- unique(agreement_sig$cohort)


# pooling results from all cohorts ----

# consider only genes where medium expression range is associated with survival
agreement_sig_range <- agreement_sig %>% 
  filter(range %in% range_of_interest)

## build contingency table

# genes where survival change and proliferation too
genes_prol_change_surv_change <- agreement_sig_range$gene[
  agreement_sig_range$agreement_surv_prol %in% c("agree", "disagree")]

# genes where survival does not change, but proliferation does
genes_prol_change_surv_not <- agreement_non_sig$gene[
  agreement_non_sig$proliferation_change == "yes"]

# genes where survival changes but proliferation does not
genes_prol_not_surv_change <- agreement_sig_range$gene[
  agreement_sig_range$agreement_surv_prol == "Ki67_proliferation_not_significant"]

# genes where neither survival nor proliferation changes
genes_prol_not_surv_not <- agreement_non_sig$gene[
  agreement_non_sig$proliferation_change == "no"]

## contigency table
m <- matrix(
  c(length(genes_prol_change_surv_change), 
    length(genes_prol_change_surv_not), 
    length(genes_prol_not_surv_change), 
    length(genes_prol_not_surv_not)), 
  nrow = 2, 
  ncol = 2, 
  byrow = T, 
  dimnames = list(
    c("proliferation_change", "proliferation_not_change"),
    c("surv_change", "surv_not_change")))

## Fisher test
fisher <- fisher.test(m)
pval <- fisher$p.value
estimate <- fisher$estimate

# result
res <- c(
  "pvalue" = pval, 
  estimate)


# medium - cohort-wise -----

# loop over cohorts and run Fisher test
fisher_cohort_l <- lapply(seq(length(cohorts)), function(i) {
  
  cat("Process", cohorts[i], "cohort\n", 
      "--------------------", 
      "\n\n")
  
  # wrangle ----
  
  # subset agreement_sig and agreement_non_sig 
  # with samples of that cohort 
  agreement_sig_cohort <- agreement_sig %>% 
    filter(cohort == cohorts[i])
  
  agreement_non_sig_cohort <- agreement_non_sig %>% 
    filter(cohort == cohorts[i])
  
  # subset agreement_sig_cohort picking only genes where range_of_interest
  # is associated with extreme prognosis
  agreement_sig_cohort_range <- agreement_sig_cohort %>% 
    filter(range %in% range_of_interest)
  
  
  # Fisher's ----
  
  # genes where survival change and proliferation too
  genes_prol_change_surv_change <- agreement_sig_cohort_range$gene[
    agreement_sig_cohort_range$agreement_surv_prol %in% c("agree", "disagree")]
  
  # non significant genes where proliferation changes
  genes_prol_change_surv_not <- agreement_non_sig_cohort$gene[
    agreement_non_sig_cohort$proliferation_change == "yes"]
  
  # genes where survival changes but proliferation does not
  genes_prol_not_surv_change <- agreement_sig_cohort_range$gene[
    agreement_sig_cohort_range$agreement_surv_prol == "Ki67_proliferation_not_significant"]
  
  # non significant genes where proliferation does not changes
  genes_prol_not_surv_not <- agreement_non_sig_cohort$gene[
    agreement_non_sig_cohort$proliferation_change == "no"]
  
  ## contigency table
  m <- matrix(
    c(length(genes_prol_change_surv_change), 
      length(genes_prol_change_surv_not), 
      length(genes_prol_not_surv_change), 
      length(genes_prol_not_surv_not)), 
    nrow = 2, 
    ncol = 2, 
    byrow = T, 
    dimnames = list(
      c("proliferation_change", "proliferation_not_change"),
      c("surv_change", "surv_not_change")))
  
  ## Fisher test
  fisher <- fisher.test(m)
  pval <- fisher$p.value
  estimate <- fisher$estimate
  
  # result
  res <- c(
    "cohort" = cohorts[i],
    "pvalue" = pval, 
    estimate)
  
  return(res)
  
})

# rbind elements of list
fisher_cohort_df <- as.data.frame(do.call(rbind, fisher_cohort_l))

# save fisher results
write.table(
  fisher_cohort_df, 
  file = file.path(
    out_dir, 
    paste0(
      "association_cohort-wise_surv-prol_", 
      paste0(range_of_interest, collapse = "-"), 
      ".tsv")
  ), 
  sep = "\t", 
  row.names = F, 
  col.names = T)
