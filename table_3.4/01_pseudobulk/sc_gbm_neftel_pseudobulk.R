# Process single-cell brain dataset from Neftel et al. 2019
# run this script in new workstation (more powerful than old one)

# packages ----
library(Seurat)
library(dplyr)


# cell ----

# expression sparse matrix
exp <- ReadMtx(
  mtx = "../analyses/single-cell/input/Neftel2019/Exp_data_TPM.mtx",
  cells = "../analyses/single-cell/input/Neftel2019/Cells.csv", 
  features = "../analyses/single-cell/input/Neftel2019/Genes.txt", 
  cell.column = 1, 
  cell.sep = ",",
  skip.cell = 1,
  feature.column = 1)

# metadata
meta <- read.csv(
  "../analyses/single-cell/input/Neftel2019/Meta-data_Neftel2019_Brain.csv")

# annotation
ann <- read.csv(
  "../analyses/single-cell/input/Neftel2019/Cells.csv")


# metadata ----

# consider only smart samples
meta_smart <- meta[meta$technology == "SmartSeq2", ]

# primary
meta_smart_primary <- meta_smart[meta_smart$sample_primary_met == "primary", ]

# type of cancer
table(meta_smart_primary$cancer_type)

# number of patients
length(table(meta_smart_primary$patient))


# annotation ----

# type of cells
table(ann$cell_type)

# possible cell types
cell_type <- names(table(ann$cell_type))


# wrangling ----

# remove cells where all genes == 0
samples_to_keep <- colSums(exp, na.rm = T) != 0

# patients to keep
patient_to_keep <- ann$sample[samples_to_keep]
table(patient_to_keep)

# remove genes where all samples == 0
genes_to_keep <- rowSums(exp, na.rm = T) != 0

# subset exp
exp_sub <- exp[genes_to_keep, samples_to_keep]

rm(exp)
gc()

# transpose
exp_t <- as.data.frame(t(exp_sub))
rm(exp_sub)
gc()

# add patient column
exp_t <- exp_t %>% 
  mutate(
    "patient" = patient_to_keep, 
    .before = colnames(exp_t)[1])

# group by patient and normalize considering only non 0 entries
exp_pseudo <- exp_t %>% 
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
  file = "/home/ieo5059/mount/local/projects/survival/output/analyses/kaplan-meier/single-cell/output/gbm_neftel/exp_pseudobulk_all-cells.rds")