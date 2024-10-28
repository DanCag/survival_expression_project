# import the processed expression of single-cell dataset (not separated by cell)
# and classify gene expression distributions into bimodal, normal-like and
# positively-skewed genes


# packages ----
library(parallel)
# library(dplyr)
source("./functions/classification_f.R")


# variables ----

# cohort
cancer_type <- "hnsc_puram"

# output dir
output_dir <- file.path(
  "../analyses/single-cell/output", 
  cancer_type)


# data ----

# single-cell

# sc dataset
df <- readRDS(
  file.path(
    "../analyses/single-cell/output", 
    cancer_type, 
    "exp_pseudobulk_all-cells.rds")
)


# loop ----

# loop across genes
classification_l <- mclapply(colnames(df)[-1], function(g) {
  
  cat("Process", g, "\n")
  
  # classify gene expression
  cl <- classification_f(gene = g, exp = df)
  
}, mc.cores = 4)

# bind rows
classification_df <- as.data.frame(do.call(rbind, classification_l))

# save classification_df_l
saveRDS(
  classification_df, 
  file = file.path(output_dir, "classification_pool-cells.rds")
)
