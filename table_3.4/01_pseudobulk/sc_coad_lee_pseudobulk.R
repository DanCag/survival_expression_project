# process single-cell coad dataset from Lee et al., 2020 (GSE132465)
# run this script in the new workstation


# packages ----
library(dplyr)


# data ----

# expression
exp_df <- read.delim(
  gzfile("../analyses/single-cell/input/GSE132465/GSE132465_GEO_processed_CRC_10X_natural_log_TPM_matrix.txt.gz"))

# annotation
ann <- read.delim(
  gzfile("../analyses/single-cell/input/GSE132465/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"))


# wrangling ----

# transform ln(TPM+1) into TPM
str(exp)
tpm <- exp(exp_df[, -1])-1

# remove normal cells
ann_tumor <- ann[ann$Class == "Tumor", ]

# number of patients
length(table(ann_tumor$Patient))

# keep only tumor cells
cells_to_keep <- ann_tumor$Index

# genes 
genes <- exp_df$Index

# transform "." into "-" in colnames of exp
colnames(tpm) <- sub("\\.", "-", colnames(tpm))

# subset expression keeping only cells_to_keep
tpm_to_keep <- tpm[, colnames(tpm) %in% cells_to_keep]
tpm_to_keep[1:3, 1:4]
rm(exp_df)

# transpose matrix
tpm_t <- t(tpm_to_keep)
colnames(tpm_t) <- genes
tpm_t[1:5, 1:5]

# cell types
cell_type <- names(table(ann_tumor$Cell_type))


# wrangling ----

# cells are in the same order in annotation and tpm_t?
all(rownames(tpm_t) == ann_tumor$Index)

# samples to keep
samples_to_keep <- rowSums(tpm_t, na.rm = T) != 0

# patients to keep
patient_to_keep <- ann_tumor$Patient[samples_to_keep]
table(patient_to_keep)

# remove genes where all samples == 0
genes_to_keep <- colSums(tpm_t, na.rm = T) != 0

# subset exp
tpm_t_sub <- as.data.frame(tpm_t[samples_to_keep, genes_to_keep])
rm(tpm_t)
gc()

# add patient column
tpm_t_sub <- tpm_t_sub %>% 
  mutate(
    "patient" = patient_to_keep, 
    .before = colnames(tpm_t_sub)[1])


# group by patient and normalize considering only non 0 entries
exp_pseudo <- tpm_t_sub %>% 
  group_by(patient) %>% 
  summarise_all(~mean(.[. > 0], na.rm = T))

# proportion of NaN
prop_nan <- sapply(
  seq(ncol(exp_pseudo))[-1], 
  function(j) {mean(is.nan(exp_pseudo[[j]]))})

# assign gene names to prop_nan
names(prop_nan) <- colnames(exp_pseudo)[-1]

# keep only genes with less than 50 % NaN
genes_to_keep <- names(prop_nan)[prop_nan < 0.5]

# # genes with less than 50 % dropout
# genes_to_keep <- colMeans(exp_patient == 0) < 0.5

# subset
exp_pseudo <- exp_pseudo[, c("patient", genes_to_keep)]

# save list
saveRDS(
  object = exp_pseudo, 
  file = "../analyses/single-cell/output/coad/exp_pseudobulk_all-cells.rds")
