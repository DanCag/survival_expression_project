# check if bimodal distributions in sc dataset with all cells pooled
# are more likely to be associated with the expression of 2 different cell types
# compared to other shapes

# packages ----
library(dplyr)


# parameters ----

# cancer type
cancer_type <- "hnsc_puram"

# threshold
threshold <- 0.05

# correction
correction <- "p_adjust_bh"


# data ----

# list of differentially expressed genes between cell types in scRNAseq
diff_exp_df_l <- readRDS(
  file.path(
    "../analyses/single-cell/output", 
    cancer_type, 
    "diff_exp_df_mean-no-zero_l.rds"))

# classification of gene expression distribution in sc dataset
classification_sc <- readRDS(
  file.path(
    "../analyses/single-cell/output", 
    cancer_type, 
    "classification_pool-cells.rds"))


# wrangling ----

# genes for which we compared expression between different cell types
genes <- diff_exp_df_l[[1]]$gene

# subset common picking only bimodal genes
classification_b <- classification_sc$gene[
  classification_sc$label == "bimodal"]

# subset common picking only non bimodal
classification_nb <- classification_sc$gene[
  classification_sc$label != "bimodal"]


# proportion of differential expressed genes ----
proportion_diff_exp <- sapply(names(diff_exp_df_l), function(x) {
  
  cat("Process", x, "\n")
  
  # dataframe with differential expression
  df <- diff_exp_df_l[[x]]
  
  res <- mean(df[[correction]] < threshold)
  
  return(res)
  
})


# association ----

# loop over cell types and run Fisher test
association_l <- lapply(names(diff_exp_df_l), function(x) {
  
  cat("Process", x, "\n")
  
  # dataframe with differential expression
  df <- diff_exp_df_l[[x]]
  
  ## bimodal diff
  bimodal_diff <- intersect(
    classification_b, 
    df$gene[df[[correction]] < threshold])
  
  # bimodal non sig
  bimodal_non_diff <- intersect(
    classification_b, 
    df$gene[df[[correction]] >= threshold])
  
  ## non bimodal diff
  non_bimodal_diff <- intersect(
    classification_nb, 
    df$gene[df[[correction]] < threshold])
  
  ## non bimodal non sig
  non_bimodal_no_diff <- intersect(
    classification_nb, 
    df$gene[df[[correction]] >= threshold])
  
  # contigency table
  m <- matrix(
    c(length(bimodal_diff), 
      length(bimodal_non_diff), 
      length(non_bimodal_diff), 
      length(non_bimodal_no_diff)), 
    nrow = 2, 
    ncol = 2, 
    byrow = T, 
    dimnames = list(
      c("bimodal", "not_bimodal"),
      c("diff_expressed", "not_diff_expressed")))
  
  fisher <- fisher.test(m)
  pval <- fisher$p.value
  estimate <- fisher$estimate
  
  # combine pvalue and estimate
  res <- c(
    "cell_type" = x,
    "pvalue" = pval,
    estimate)
  
  return(res)
  
})

# bind elements of list into dataframe
association_df <- as.data.frame(do.call(rbind, association_l))

# add column with proportion of differential expressed genes
association_df <- association_df %>% 
  mutate("prop_diff_exp" = proportion_diff_exp, .after = cell_type)

# output directory
out_dir <- file.path(
  "../output/analyses/kaplan-meier/association/single-cell/mean-no-zero/compare-within-sc", 
  cancer_type,
  correction,
  paste0(
    "qvalue-",
    as.character(threshold)))

# create directory(ies) recursively if not there
dir.create(
  out_dir,
  showWarnings = FALSE ,
  recursive = T)

# save
write.table(
  association_df, 
  file = file.path(
    out_dir,
    paste0(
      "bimodal_diff-expression_", 
      correction,
      "_mean-no-zero_",
      "qvalue-",
      as.character(threshold), 
      ".tsv")), 
  sep = "\t", 
  row.names = F, 
  col.names = T)
