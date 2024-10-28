# Extract survival tables from mysql database
# Process survival tables
# Merge survival tables with expression matrices

# Run this script on workstation

# import packages ----
library(RMySQL)
library(tidyverse)


# database ----

# set up connection with database
# you need to have access to database
mydb = dbConnect(
  MySQL(), 
  user='username', 
  password='password', 
  dbname='tcga', 
  host='146.247.16.202')

# show available tables
dbListTables(mydb)


# input variables ----

# survival analysis 
surv_analysis <- "DSS"

# imputation info
imp_info <- "no-imp"


# fetch data ----

# query table with columns' names of CBIO
rs <- dbSendQuery(mydb, 
                  "SELECT * FROM tcga.CBIO_clinical_NAMES;")

# fetch data from database
cbioportal_clinical_columns_names <- fetch(rs, n = -1)

# query cbioportal table
query <- dbSendQuery(mydb, 
                     "SELECT * FROM tcga.CBIO_clinical;")

# fetch data from database
cbioportal_clinical <- fetch(query, n = -1)


# Data wrangling ----

# according to the type of survival analysis, 
# there will be different columns to be selected
if (surv_analysis == "DSS") {
  
  # vector of surv columns to select
  surv_cols <- c(
    "Disease-specific Survival status",
    "Months of disease-specific survival")
  
} else if (surv_analysis == "DFS") {
  
  # vector of surv columns to select
  surv_cols <- c(
    "Disease Free Status",
    "Disease Free (Months)")
  
} else if (surv_analysis == "PFS") {
  
  # vector of surv columns to select
  surv_cols <- c(
    "Progression Free Status",
    "Progress Free Survival (Months)")
  
}

# vector of columns to select
cols_to_select <- c(
  "TCGA PanCanAtlas Cancer Type Acronym", 
  "Patient ID",
  "Sample ID",
  surv_cols)

# subset cbioportal clinical by picking 
# only primary samples and selecting specific columns
cbioportal_clinical_sub <- cbioportal_clinical[
  cbioportal_clinical$`Sample Type` == "Primary", 
  cols_to_select]

# remove NAs
cbioportal_clinical_sub_no_NA <- na.omit(cbioportal_clinical_sub)

# split status column at ":"
if (surv_analysis == "DSS") {
  
  cbioportal_clinical_sub_split <- cbind(
    cbioportal_clinical_sub_no_NA[, c(1, 2, 3)], 
    do.call('rbind', 
            strsplit(
              cbioportal_clinical_sub_no_NA[["Disease-specific Survival status"]], ':', fixed=TRUE)),
    cbioportal_clinical_sub_no_NA[, 5]) 
  
} else if (surv_analysis == "DFS") {
  
  cbioportal_clinical_sub_split <- cbind(
    cbioportal_clinical_sub_no_NA[, c(1, 2, 3)], 
    do.call('rbind', 
            strsplit(
              cbioportal_clinical_sub_no_NA[["Disease Free Status"]], ':', fixed=TRUE)),
    cbioportal_clinical_sub_no_NA[, 5]) 
  
} else if (surv_analysis == "PFS") {
  
  cbioportal_clinical_sub_split <- cbind(
    cbioportal_clinical_sub_no_NA[, c(1, 2, 3)], 
    do.call('rbind', 
            strsplit(
              cbioportal_clinical_sub_no_NA[["Progression Free Status"]], ':', fixed=TRUE)),
    cbioportal_clinical_sub_no_NA[, 5]) 
}

# rename col 4 and 5
colnames(cbioportal_clinical_sub_split) <- c(
  colnames(cbioportal_clinical_sub_split)[1:3], 
  paste("event", "int", sep = "_"), 
  paste("event", "chr", sep = "_"), 
  "months")

# make column 4 an integers column
class(cbioportal_clinical_sub_split[[4]]) <- "integer"

# add column with years
cbioportal_clinical_sub_split$years <- cbioportal_clinical_sub_split$months/12

# all TCGA cohorts (33)
cohorts <- names(
  table(
    cbioportal_clinical_sub_split$`TCGA PanCanAtlas Cancer Type Acronym`)
)

# subset TCGA cohorts taking the ones in common between 
# cbioportal_clinical_sub_split and cohorts
cohorts <- names(table(cbioportal_clinical_sub_split$`TCGA PanCanAtlas Cancer Type Acronym`))


# loop ----

# loop over cohorts and extract survival info
l <- lapply(cohorts, function(x) {
  
  cat("Process", x, "cohort","\n", 
      "-----------------", 
      "\n\n")
  
  # pick only patients from desired cancer cohort
  # and desired columns
  cbioportal_cohort <- cbioportal_clinical_sub_split %>% 
    filter(`TCGA PanCanAtlas Cancer Type Acronym` == x) 
  
  # gene expression matrix path 
  # Data are in the cluster. You need to have the cluster directory mounted
  # in this case I have mounted the cluster on the workstation
  gene_exp_path <- file.path(
    "/home/ieo5059/mount/cl/data/MS/tcga_open",
    "all_tcga_projects_fpkm_and_tpm_from_gdc",
    "gdc_tcga_all_projects_tpm_tables_with_mitochondrial_genes", 
    paste0("TCGA-", 
           x, 
           "_tpm_df.txt"))
  
  
  # import datasets ----
  
  cat("Import gene expression dataset", 
      "\n\n")
  
  # gene expression
  exp <- read.delim(gene_exp_path, 
                    check.names = F)
  
  
  # Data wrangling ----
  
  # vector with trimmed ensembl_id (no version)
  ensembl_id_trim <- sub("\\..*", "", exp$ensembl_id)
  
  # transpose dataframe
  exp_t <- as.data.frame(t(exp[, -1]))
  colnames(exp_t) <- ensembl_id_trim
  
  # get sample id
  exp_t$sample_id <- substr(colnames(exp)[2:ncol(exp)], 1, 15)
  
  # if there is any duplicated samples
  if (any(duplicated(exp_t$sample_id))) {
    
    cat("Collapse same samples into one by taking the median", 
        "\n\n")
    
    # take the median expression value for each gene
    exp_t <- exp_t %>% 
      group_by(sample_id) %>% 
      summarise_all(median)
    
  }
  
  # merge survival table with expression matrix
  surv_exp <- merge(cbioportal_cohort, 
                    exp_t,
                    by.x = "Sample ID", 
                    by.y = "sample_id"
  )
  
  # save ----
  # save table with survival time and expression 
  
  cat("Save survival-expression table", 
      "\n\n")
  
  # output directory (mounted local machine on the workstation)
  output_dir <- file.path(
    "/home/ieo5059/mount/local/projects/survival/output/survival-tables/process",
    paste("surv-time_events-censored", imp_info, sep ="_"),
    "survival-exp",
    "continuos_exp",
    x, 
    surv_analysis)
  
  # output file
  outfile <- file.path(
    output_dir, 
    paste0("surv-time_exp_events-censored_", 
           imp_info, 
           "_", 
           surv_analysis,
           "_",
           x,
           ".rds"))
  
  
  # create directory(ies) recursively if not there
  dir.create(output_dir, 
             showWarnings = FALSE ,
             recursive = T)
  
  saveRDS(object = surv_exp, 
          file = outfile)
  
  cat("\n")
})
