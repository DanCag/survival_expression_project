# How does distribution shape of middle extreme genes spotted in TCGA cohorts
# change in healthy samples (GTEx)?

# example of contingency table
#                       # medium-extreme genes    # non significant genes
# shape changes
# shape does not change


# packages ----
library(dplyr)


# input parameters ----

# cancer cohorts
cohorts <- readRDS("../data/cohorts.rds")

# survival analysis
surv_analysis <- "DSS"

# significance threshold for multiple testing correction
qt <- 0.2

# significance threshold for log-rank test pairwise comparisons
kmpt <- 0.05

# list to match cohorts with GTEx tissues
cohort_tissue_l <- list(
  "adrenal_gland" = "ACC", 
  "bladder" = "BLCA", 
  "bone_marrow" = c("LAML", "LCML"), # useless because no expression for LAML
  "brain" = "GBM", 
  "breast" = "BRCA", 
  "cervix_uteri" = "CESC", 
  "colon" = "COAD", 
  "esophagus" = "ESCA", 
  "kidney" = c("KICH", "KIRC", "KIRP"), 
  "liver" = c("LIHC"),
  "lung" = c("LUAD", "LUSC"),
  "muscle" = "SARC", 
  "ovary" = "OV", 
  "pancreas" = "PAAD", 
  "prostate" = "PRAD", 
  "skin" = "SKCM", 
  "stomach" = "STAD", 
  "testis" = "TGCT", 
  "thyroid" = "THCA", 
  "uterus" = c("UCS", "UCEC"))

# subset cohorts picking only those with a match in GTEx
cohorts_matched <- cohorts[cohorts %in% unlist(cohort_tissue_l)] 


# input/output ----

# classification - tumor
classification_tumor_dir <- "../analyses/classification/TCGA"

# classification - normal
classification_normal_dir <- "../analyses/classification/GTEx"

# count directory (genes that are significantly associated with survival)
count_dir <- file.path(
  "../analyses/counts/new/shape-wise",
  paste0(as.character(qt), "-qvalue-threshold"),
  paste0(as.character(kmpt), "-km-pairwise-threshold"))

# output 
out_dir <- "../analyses/association/tumor_normal"


# loop ----
fisher_l <- lapply(seq(cohorts_matched), function(i) {
  
  cat("Process", cohorts_matched[i], "\n")
  
  # classification - tumor - cohort
  classification_tumor_cohort_path <- file.path(
    classification_tumor_dir, 
    cohorts_matched[i], 
    paste0(
      "genes_classification_", 
      cohorts_matched[i], 
      ".rds"))
  
  # classification - normal - tissue
  classification_normal_tissue_path <- file.path(
    classification_normal_dir, 
    names(cohort_tissue_l[grep(cohorts_matched[i], cohort_tissue_l)]), 
    paste0(
      "genes_classification_", 
      names(cohort_tissue_l[grep(cohorts_matched[i], cohort_tissue_l)]), 
      ".rds"))
  
  # genes significantly associated with extreme prognosis  
  count_path <- file.path(
    count_dir, 
    cohorts_matched[i], 
    surv_analysis,
    paste0(
      "counts_", 
      cohorts_matched[i], 
      "_", 
      surv_analysis, 
      ".rds"))
  
  # check if counts file exists
  if ( file.exists(classification_tumor_cohort_path) &
       file.exists(classification_normal_tissue_path) &
       file.exists(count_path) ) {
    
    # data ----
    
    # classification in tumor cohort
    classification_tumor <- readRDS(classification_tumor_cohort_path)
    
    # classification in normal tissue
    classification_normal <- readRDS(classification_normal_tissue_path)
    
    # count 
    count <- readRDS(count_path)
    
    
    # wrangle classification ----
    
    # rename columns of tumor samples
    colnames(classification_tumor)[2:3] <- c(
      "classification_tumor", 
      "pvalue_tumor")
    
    # rename columns of normal samples classification
    colnames(classification_normal)[2:3] <- c(
      "classification_normal", 
      "pvalue_normal")
    
    # merge classification tumor and classification normal
    # we cannot define the gene expression distribution shape
    classification_tumor_normal <- merge(
      classification_tumor, 
      classification_normal, 
      by = "gene") 
    
    
    # wrangle count ----
    
    # do not consider to_disentangle hits
    count_sub <- count %>% 
      filter(
        count_low != "to_disentangle", 
        count_medium != "to_disentangle",
        count_high != "to_disentangle")
    
    # transform count columns into numeric
    count_sub <- count_sub %>% 
      mutate_at(c(
        "count_low",
        "count_medium",
        "count_high"), as.numeric)
    
    # keep genes where in medium we see an extreme prognosis 
    count_to_keep <- count_sub[
      count_sub$count_medium %in% c(-1, 1), ]
    
    # subset count keeping only those genes that are present both
    # in classification tumor and classification normal
    count_to_keep <- count_to_keep[
      count_to_keep$gene %in% classification_tumor_normal$gene, ]
    
    # remove genes associated with extreme prognosis in either low and high
    # from classification_tumor_normal in order to run the comparison only against  
    # non significant genes
    to_remove_genes <- count_sub$gene[
      count_sub$count_low %in% c(-1, 1) |
        count_sub$count_high %in% c(-1, 1)]
    
    # subset classification_tumor_normal
    classification_tumor_normal <- classification_tumor_normal[
      !classification_tumor_normal$gene %in% to_remove_genes, ]
    
    
    # fisher ----
    
    # genes where shape changes between tumor and normal 
    # and that are associated with extreme prognosis in medium
    shape_change_associated <- intersect(
      classification_tumor_normal$gene[
        classification_tumor_normal$classification_tumor != 
          classification_tumor_normal$classification_normal], 
      count_to_keep$gene)
    
    # genes where shape changes between tumor and normal
    # and that are not significant
    shape_change_not_associated <- setdiff(
      classification_tumor_normal$gene[
        classification_tumor_normal$classification_tumor != 
          classification_tumor_normal$classification_normal], 
      count_to_keep$gene)
    
    # genes where shape does not change between tumor and normal
    # and that are associated with extreme prognosis in medium
    shape_not_change_associated <- intersect(
      classification_tumor_normal$gene[
        classification_tumor_normal$classification_tumor == 
          classification_tumor_normal$classification_normal], 
      count_to_keep$gene)
    
    # genes where shape does not change between tumor and normal
    # and that are not significant
    shape_not_change_not_associated <- setdiff(
      classification_tumor_normal$gene[
        classification_tumor_normal$classification_tumor == 
          classification_tumor_normal$classification_normal], 
      count_to_keep$gene)
    
    # contingency table
    m <- matrix(
      c(length(shape_change_associated), 
        length(shape_change_not_associated), 
        length(shape_not_change_associated), 
        length(shape_not_change_not_associated)), 
      nrow = 2, 
      ncol = 2, 
      byrow = T, 
      dimnames = list(
        c("shape_changes", "shape_not_change"),
        c("associated", "not_associated")))
    
    fisher <- fisher.test(m)
    pval <- fisher$p.value
    estimate <- fisher$estimate
    
    # chisquare
    # chisq_res <- chisq.test(m)
    
    # combine pvalue and estimate
    res <- c(
      cohorts_matched[i],
      pval,
      estimate)
    
  } else {
    
    cat("No files for ", cohorts_matched[i], "\n")
    
    res <- c(
      cohorts_matched[i],
      rep(NA, 2))
    
  }
  
  names(res) <- c(
    "cohort", 
    "pvalue", 
    "odds_ratio")
  
  return(res)
  
})

# combine elements in info_l
fisher_df <- bind_rows(fisher_l)

# transform columns into numeric
fisher_df <- fisher_df %>% 
  mutate_at(c(2:3), as.numeric)

# save table with associations
write.table(
  fisher_df, 
  file = file.path(
    out_dir, 
    "association_shape-change_middle-extreme-prognosis.tsv"), 
  sep = "\t", 
  row.names = F, 
  col.names = T)
