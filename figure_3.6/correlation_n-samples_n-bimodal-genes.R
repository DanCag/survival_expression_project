# correlation between number of samples and fraction of bimodal genes
# either in normal (GTEx) or tumor (TCGA) samples


# packages ----
library(dplyr)
library(ggplot2)
library(ggpubr)


# input parameters ----

# cancer cohorts
cohorts <- readRDS("../data/cohorts.rds")

# survival analysis
surv_analysis <- "DSS"

# significance threshold for multiple testing correction
qt <- 0.2

# significance threshold for log-rank test pairwise comparisons
kmpt <- 0.05

# interested in tumor or normal?
interest <- "tumor"

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

# depending if we deal with normal or tumor samples, cohorts change
if (interest == "normal") {
  
  # subset cohorts picking only those with a match in GTEx
  cohorts <- cohorts[cohorts %in% unlist(cohort_tissue_l)]
  
} else {
  
  cohorts <- cohorts
}


# input/output ----

# input
if (interest == "normal") {
  
  # classification dir
  classification_dir <- "../analyses/classification/GTEx"
  
  # number of samples file 
  n_samples <- read.delim(
    "../analyses/patients/new/tcga_gtex_n-samples.tsv", 
    stringsAsFactors = F, 
    check.names = F)[, c(1:2)]
  
} else {
  
  # classification dir
  classification_dir <- "../analyses/classification/TCGA"
  
  # number of samples file 
  n_samples <- read.delim(
    "../analyses/patients/new/tcga_normal-tumor_n-samples.tsv", 
    stringsAsFactors = F, 
    check.names = F)[, c(1, 3)]

}


# loop ----

# loop across cohorts
n_freq_l <- lapply(seq(cohorts), function(i) {
  
  cat("Process", cohorts[i], "\n")
  
  if (interest == "normal") {
    
    # classification - tissue
    classification_tissue_path <- file.path(
      classification_dir, 
      names(cohort_tissue_l[grep(cohorts[i], cohort_tissue_l)]), 
      paste0(
        "genes_classification_", 
        names(cohort_tissue_l[grep(cohorts[i], cohort_tissue_l)]), 
        ".rds"))
    
  } else {
    
    # classification - tissue
    classification_tissue_path <- file.path(
      classification_dir, 
      cohorts[i], 
      paste0(
        "genes_classification_", 
        cohorts[i], 
        ".rds"))
    
  }
  
  # check if counts file exists
  if ( file.exists(classification_tissue_path) ) {
    
    
    # data ----
    
    # classification 
    classification <- readRDS(classification_tissue_path)
    
    if (interest == "normal") {
      
      # number of samples
      n <- n_samples[[paste0("n_samples_", interest)]][
        grep(names(cohort_tissue_l[grep(cohorts[i], cohort_tissue_l)]), 
             n_samples$tissue)]
      
    } else {
      
      # number of samples
      n <- n_samples[[paste0("n_samples_", interest)]][
        grep(cohorts[i], n_samples$cohort)]
    }
    
    # if classification table is not empty
    if ( nrow(classification) != 0 ) {
      
      # frequency of bimodal genes
      freq <- mean(classification$label == "bimodal")
      
    } else {
      
      # we cannot compute frequency
      freq <- NA
    }
    
    # result
    res <- c(
      "cohort" = cohorts[i], 
      "n_samples" = n, 
      "frequency" = freq)
    
  } else {
    
    cat("No files for ", cohorts[i], "\n")
    
    res <- c(
      "cohort" = cohorts[i],
      setNames(
        object = rep(NA, 2), 
        nm = c("n_samples", "frequency"))
    )
    
  }
  
  return(res)
  
})

# combine elements in info_l
n_freq_df <- bind_rows(n_freq_l)

# transform columns into numeric
n_freq_df <- n_freq_df %>% 
  mutate_at(c(2:3), as.numeric)


# correlation ----

# correlation between number of samples and frequency of bimodal genes
cor.test(
  n_freq_df$n_samples,
  n_freq_df$frequency,
  method = "spearman")

# scatterplot
ggplot(n_freq_df, aes(x = n_samples, y = frequency)) + 
  geom_point(size = 3) + 
  ggrepel::geom_text_repel(
    label = n_freq_df$cohort) +
  stat_cor(method = "spearman") +
  theme_bw()