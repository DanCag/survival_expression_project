# compare number of bimodal middle extreme genes I get when
# - using tertiles split
# - using shape-wise split

# packages ----
library(dplyr)
library(ggplot2)


# parameters ----

# survival analysis
surv_analysis <- "DSS"

# significance threshold for multiple testing correction 
qt <- 0.2

# significance threshold for pairwise KM comparison
kt <- 0.05

# method for multiple test correction
method <- "BH"


# input ----

# classification directory
classification_dir <- "../analyses/classification/TCGA"

# count directory tertiles
count_tertiles_dir <- file.path(
  "../analyses/counts/new/tertiles",
  paste0(
    as.character(qt),
    "-qvalue-threshold"), 
  paste0(
    as.character(kt), 
    "-km-pairwise-threshold"))

# count directory shape-wise
count_shape_wise_dir <- file.path(
  "../analyses/counts/new/shape-wise",
  paste0(
    as.character(qt),
    "-qvalue-threshold"), 
  paste0(
    as.character(kt), 
    "-km-pairwise-threshold"))

# cancer cohorts
cohorts <- readRDS("../data/cohorts.rds")


# loop ----

# loop over cohorts
ratio_l <- lapply(seq(length(cohorts)), function(i) {
  
  cat("Process", cohorts[i], "cohort\n", 
      "--------------------", 
      "\n\n")
  
  # classification file
  classification_path <-  file.path(
    classification_dir, 
    cohorts[i], 
    paste0(
      "genes_classification_", 
      cohorts[i], 
      ".rds"))
  
  # count tertiles path
  count_tertiles_path <- file.path(
    count_tertiles_dir, 
    cohorts[i], 
    surv_analysis, 
    paste0(
      "counts_", 
      cohorts[i],
      "_", 
      surv_analysis,
      ".rds"))
  
  # count shape-wise path
  count_shape_wise_path <- file.path(
    count_shape_wise_dir, 
    cohorts[i], 
    surv_analysis, 
    paste0(
      "counts_", 
      cohorts[i],
      "_", 
      surv_analysis,
      ".rds"))  
  
  
  # check if files exist
  if ( file.exists(classification_path) &
       file.exists(count_tertiles_path) & 
       file.exists(count_shape_wise_path) ) {
    
    # data ----
    
    # import dataframes
    classification <- readRDS(classification_path)
    
    # import count tertiles
    count_tertiles <- readRDS(count_tertiles_path)
    
    # import count shape-wise
    count_shape_wise <- readRDS(count_shape_wise_path)
    
    
    # wrangling ----
    
    # bimodal genes
    bimodal_genes <- classification$gene[classification$label == "bimodal"]
    
    # tertiles - keep only bimodal middle extreme genes
    tertiles_bimodal_me <- count_tertiles$gene[
      count_tertiles$gene %in% bimodal_genes & 
        count_tertiles$count_medium %in% c(-1, 1)]
    
    # shape-wise - keep only bimodal middle extreme genes
    shape_wise_bimodal_me <- count_shape_wise %>% 
      filter(label == "bimodal", count_medium %in% c(-1, 1)) %>% 
      pull(gene)
    
    # ratio tertiles-bimodal-me/shape-wise-bimodal-me
    ratio <- length(tertiles_bimodal_me)/length(shape_wise_bimodal_me)
    
    # result
    res <- c(
      "cohort" = cohorts[i],
      "n_tertiles_bimodal_me" = length(tertiles_bimodal_me),
      "n_shape_wise_bimodal_me" = length(shape_wise_bimodal_me),
      "ratio" = ratio)
    
  } else {
    
    cat("No files for ", cohorts[i],"\n")
    res <- NULL
  }
  
  return(res)
})

# rbind elements of list
ratio_df <- as.data.frame(do.call(rbind, ratio_l))  

# save table
write.table(
  ratio_df, 
  file = "../analyses/tertile_shape-wise_count-comparison/me_count-comparison_bimodal.tsv",
  sep = "\t", 
  col.names = T, 
  row.names = F)


# barplot ----

# remove cohorts where in both case we observe 0 genes
# remove cohorts where there is no change between the two splits
ratio_df_sub <- ratio_df[
  (ratio_df$n_tertiles_bimodal_me != 0 & ratio_df$n_shape_wise_bimodal_me != 0) &
    (ratio_df$n_tertiles_bimodal_me != ratio_df$n_shape_wise_bimodal_me), -ncol(ratio_df)]

# rename columns
colnames(ratio_df_sub)[2:3] <- c("tertile", "shape_wise") 

# gather
ratio_df_sub_gather <- ratio_df_sub %>% 
  tidyr::gather("split", "n", -cohort)

# factor 
ratio_df_sub_gather$split <- factor(ratio_df_sub_gather$split, levels = c("tertile", "shape_wise"))

# barplot
ggplot(
  ratio_df_sub_gather, aes(x = forcats::fct_rev(cohort), y = n, fill = forcats::fct_rev(split))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  coord_flip()