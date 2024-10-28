# label distribution of healthy samples (GTEx)
# run this script in the worstation


# packages ----
library(parallel)
library(dplyr)
source("./functions/classification_f.R")


# parameters ----

# expression path
exp_path <- "../data/expression_normal_samples/GTEx/GTEx_v8_exp.rds"

# annotation path
ann_path <- "../data/expression_normal_samples/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

# tcga genes
tcga_genes_path <- "../analyses/genes_list/tcga_genes_ensembl.rds"

# number of cores
n_cores <- 20


# output ----

# directory where we will store distribution shape classification
classification_dir <- "../analyses/classification/new/GTEx"


# import data ----

# gtex expression
gtex_exp <- readRDS(exp_path)

# gtex annotation
gtex_ann <- read.delim(
  ann_path, 
  as.is = TRUE, 
  header = TRUE)

# tcga genes 
tcga_genes <- readRDS(tcga_genes_path)


# wrangling ----

# explore tissues 
table(gtex_ann$SMTS)

# vector with tissues of interest
tissues <- c(
  "Adrenal Gland", 
  "Bladder", 
  "Bone Marrow", 
  "Brain", 
  "Breast", 
  "Cervix Uteri", 
  "Colon", 
  "Esophagus", 
  "Kidney", 
  "Liver",
  "Lung",
  "Muscle", 
  "Ovary", 
  "Pancreas", 
  "Prostate", 
  "Skin", 
  "Stomach", 
  "Testis", 
  "Thyroid", 
  "Uterus")

# subset annotation
gtex_ann_sub <- gtex_ann[gtex_ann$SMTS %in% tissues, ]
table(gtex_ann_sub$SMTS)

# genes in common with tcga
common_genes <- intersect(
  tcga_genes, 
  sub("\\..*", "", gtex_exp$Name))


# loop ----

# loop over GTEx tissue and classify gene expression
gtex_info <- sapply(seq(tissues), function(i) {
  
  cat("Process", tissues[i], "\n")
  
  # samples of specific tissue
  samples_to_select <- gtex_ann$SAMPID[gtex_ann$SMTS == tissues[i]]
  
  # subset the expression picking only samples_to_select
  gtex_exp_sub <- gtex_exp[, colnames(gtex_exp) %in% c("Name", "Description", samples_to_select)]
  
  # transpose gtex_exp_sub and assign ensembl to colnames
  gtex_exp_sub_t <- t(gtex_exp_sub[, -c(1,2)])
  colnames(gtex_exp_sub_t) <- gtex_exp_sub$Name
  
  
  # loop ----
  
  # loop over genes and classify distributions
  
  # rename col names of gtex_exp_sub_t removing the version info
  colnames(gtex_exp_sub_t) <- sub("\\..*", "", colnames(gtex_exp_sub_t))
  
  # subset gtex expression
  gtex_exp_common <- as.data.frame(gtex_exp_sub_t[, common_genes])
  
  # remove genes with all 0s 
  cols_exp_selected <- colnames(gtex_exp_common)[
    apply(gtex_exp_common, 
          2, 
          function(x) sum(x) != 0)]
  
  
  # loop ----
  
  # loop over genes and label them
  
  # loop over genes and label them
  distr_classification_l <- mclapply(
    cols_exp_selected,
    function(g) {
      
      # cat("Assess", g, "gene distribution\n")
      cl <- classification_f(gene = g, exp = gtex_exp_common)
      
    }, mc.cores = n_cores)
  
  # dataframe with gene, classification and pvalue information
  classification_df <- bind_rows(distr_classification_l)
  
  # classification path
  classification_tissue_dir <- file.path(
    classification_dir,
    tolower(sub(" ", "_", tissues[i])))
  
  # create directory(ies) recursively if not there
  dir.create(
    classification_tissue_dir,
    showWarnings = FALSE ,
    recursive = T)
  
  # save table
  saveRDS(
    classification_df, 
    file = file.path(
      classification_tissue_dir,
      paste0(
        "genes_classification_",
        tolower(sub(" ", "_", tissues[i])),
        ".rds")))
  
  return(dim(gtex_exp_common)[1])
})

# make dataframe with info about the number of samples in each cohort
gtex_info_df <- data.frame(
  "cohort" = tolower(sub(" ", "_", tissues)), 
  gtex_info)

# save info
write.table(
  gtex_info_df, 
  file = "../analyses/patients/new/gtex_n-samples.tsv", 
  sep = "\t", 
  row.names = F, 
  col.names = T)