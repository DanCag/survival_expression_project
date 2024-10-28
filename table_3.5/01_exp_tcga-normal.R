# Extract TCGA normal samples expression from mysql database
# run this script in the workstation

# fix this script

# import packages ----
library(RMySQL)
library(tidyverse)


# database ----

# set up db connection
mydb = dbConnect(
  MySQL(), 
  user='username', 
  password='password', 
  dbname='tcga', 
  host='146.247.16.202')

# show available tables
dbListTables(mydb)


# input variables ----

# cancer cohorts
cohorts <- readRDS("../data/cohorts.rds")


# output ----

# output directory where I store expression of normal samples
out_exp_dir <- "../data/expression_normal_samples/tcga_normal" 


# fetch data ----

# query TCGA TPM table picking only normal (11) samples for all cohorts
# normal sample info got from "TCGA_ids" table
query_normal = dbSendQuery(
  mydb, 
  "SELECT 
  project_id, 
  sample_submitter_id, 
  ensmbl_id, 
  tpm_counts 
  FROM 
  TCGA_tpm_counts_full t1 INNER JOIN TCGA_ids t2 ON t1.`Sample ID` = t2.sample_submitter_id
  WHERE t2.sample_type_id = 11;")

# fetch data
tpm_normal_melt <- fetch(query_normal, n = -1)


# loop ----

# loop over cohorts and get normal samples
normal_info_l <- lapply(seq(cohorts), function(j) {
  
  cat("Process", cohorts[j], "cohort\n") 
  
  # extract expression of normal samples from one cohort
  normal <- tpm_normal_melt %>% 
    filter(project_id == paste0("TCGA-", cohorts[j]))
  
  # check if there are any normal sample
  if ( nrow(normal) == 0 ) {
    
    cat("No normal samples for ", cohorts[j], "\n")
    
    # result
    info <- c(
      "cohort" = cohorts[j], 
      "n_normal_samples" = 0
    )
    
  } else {
    
    # spread dataset
    normal_wide <- spread(
      data = normal,
      key = ensmbl_id, 
      tpm_counts)
    
    # output path dir
    out_exp_cohort_dir <- file.path(
      out_exp_dir, 
      cohorts[j])
    
    # create directory(ies) recursively if not there
    dir.create(
      out_exp_cohort_dir,
      showWarnings = FALSE ,
      recursive = T)
    
    # save dataframe
    saveRDS(
      object = normal_wide[, -1], 
      file = file.path(
        out_exp_cohort_dir,
        paste0(
          "tcga_normal_", 
          cohorts[j], 
          ".rds")
      )
    )
    
    # result
    info <- c(
      "cohort" = cohorts[j], 
      "n_normal_samples" = nrow(normal_wide)
    )
    
  }
  
  return(info)
  
})

# bind elements of list 
normal_info_df <- bind_rows(normal_info_l)

# save this info
write.table(
  normal_info_df, 
  file = "../data/expression_normal_samples/tcga_normal/normal_samples_info.tsv", 
  sep = "\t", 
  row.names = F)
