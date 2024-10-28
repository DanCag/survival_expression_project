# Process single-cell pdca dataset from Hwang et al. 2022
# Run this script in new workstation (more powerful than old workstation)


# packages ----
library(Seurat)
library(dplyr)


# cell ----

# import expression matrices of different samples in dataset
# samples are organized based on different type of breast cancer

# unit of measure is transcript per 10 K

# cell 1
cell1 <- ReadMtx(
  "../analyses/single-cell/input/Hwang2022/Exp_data_TP10K_1.mtx", 
  cells = "../analyses/single-cell/input/Hwang2022/Cells1.csv", 
  features = "../analyses/single-cell/input/Hwang2022/genes.txt", 
  cell.column = 1, 
  feature.column = 1,
  cell.sep = ",",
  feature.sep = " ", 
  skip.cell = 1)


# cell 2
cell2 <- ReadMtx(
  "../analyses/single-cell/input/Hwang2022/Exp_data_TP10K_2.mtx", 
  cells = "../analyses/single-cell/input/Hwang2022/Cells2.csv", 
  features = "../analyses/single-cell/input/Hwang2022/genes.txt", 
  cell.column = 1, 
  feature.column = 1,
  cell.sep = ",",
  feature.sep = " ", 
  skip.cell = 1)

# cell 3
cell3 <- ReadMtx(
  "../analyses/single-cell/input/Hwang2022/Exp_data_TP10K_3.mtx", 
  cells = "../analyses/single-cell/input/Hwang2022/Cells3.csv", 
  features = "../analyses/single-cell/input/Hwang2022/genes.txt", 
  cell.column = 1, 
  feature.column = 1,
  cell.sep = ",",
  feature.sep = " ", 
  skip.cell = 1)


# meta ----

# metadata
meta <- read.csv(
  "../analyses/single-cell/input/Hwang2022/Meta-data.csv")

# annotation ----

# import annotation tables of different samples in dataset
# samples are organized based on different type of breast cancer

# annotation1
ann1 <- read.csv(
  "../analyses/single-cell/input/Hwang2022/Cells1.csv")

# annotation2
ann2 <- read.csv(
  "../analyses/single-cell/input/Hwang2022/Cells2.csv")

# annotation3
ann3 <- read.csv(
  "../analyses/single-cell/input/Hwang2022/Cells3.csv")


# wrangling ----

# combine annotation
ann_all <- do.call(rbind, list(ann1, ann2, ann3))

# type of cells
table(ann_all$cell_type)

# combine all cells
cell_all <- do.call(cbind, list(cell1, cell2, cell3))

# remove single subsets to free up space
rm(cell1, cell2, cell3)

# check if cell names in ann_all_sub coincide with cell_sub colnames
all(ann_all$cell_name == cell_all@Dimnames[[2]])


# cell types ids ----

# possible cell types
cell_type <- names(table(ann_all$cell_type))


# wrangling ----

# remove cells where all genes == 0
samples_to_keep <- colSums(cell_all, na.rm = T) != 0

# patients to keep
patient_to_keep <- ann_all$sample[samples_to_keep]
table(patient_to_keep)

# remove genes where all samples == 0
genes_to_keep <- rowSums(cell_all, na.rm = T) != 0

# subset exp
exp_sub <- cell_all[genes_to_keep, samples_to_keep]

rm(cell_all)
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

# subset
exp_pseudo <- exp_pseudo[, c("patient", genes_to_keep)]

# save list
saveRDS(
  object = exp_pseudo, 
  file = "../analyses/single-cell/output/paad_hwang/exp_pseudobulk_all-cells.rds")