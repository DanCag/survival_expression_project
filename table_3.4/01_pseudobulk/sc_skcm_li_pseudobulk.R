# Process single-cell skcm dataset from Li et al. 2019
# Run this script in new workstation (more powerful than old workstation)


# packages ----
library(Seurat)
library(dplyr)


# cell ----

# expression sparse matrix
exp <- ReadMtx(
  mtx = "../analyses/single-cell/input/Li2019/Exp_data_UMIcounts.mtx",
  cells = "../analyses/single-cell/input/Li2019/Cells.csv", 
  features = "../analyses/single-cell/input/Li2019/Genes.txt", 
  cell.column = 1, 
  cell.sep = ",",
  skip.cell = 1,
  feature.column = 1)

# metadata
meta <- read.csv(
  "../analyses/single-cell/input/Li2019/Meta-data.csv")

# annotation
ann <- read.csv(
  "../analyses/single-cell/input/Li2019/Cells.csv")


# metadata ----

# cancer sample
table(meta$cancer_type)

# number of patients
length(table(meta$patient))


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

# transform UMI counts into counts per million
exp_cpm <- as.data.frame(t(apply(exp_t, 1, function(x)(x/sum(x))*1e06)))
rm(exp_t)
gc()

# add patient column
exp_cpm <- exp_cpm %>% 
  mutate(
    "patient" = patient_to_keep, 
    .before = colnames(exp_cpm)[1])

# group by patient and normalize considering only non 0 entries
exp_pseudo <- exp_cpm %>% 
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

# subset
exp_pseudo <- exp_pseudo[, c("patient", genes_to_keep)]

# save list
saveRDS(
  object = exp_pseudo, 
  file = "../analyses/single-cell/output/skcm_li/exp_pseudobulk_all-cells.rds")