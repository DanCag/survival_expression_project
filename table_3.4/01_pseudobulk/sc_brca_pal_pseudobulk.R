# process single-cell breast dataset from Pal et al. 2021
# run this script in new workstation (more powerful than old one)


# packages ----
library(Seurat)
library(dplyr)


# cell ----

# import expression matrices of different samples in dataset
# samples are organized based on different type of breast cancer

# cell 1
cell1 <- ReadMtx(
  "../analyses/single-cell/input/Pal2021/Exp_data_UMIcounts1.mtx", 
  cells = "../analyses/single-cell/input/Pal2021/Cells1.csv", 
  features = "../analyses/single-cell/input/Pal2021/Genes.txt", 
  cell.column = 1, 
  feature.column = 1,
  cell.sep = ",",
  feature.sep = " ", 
  skip.cell = 1)


# cell 2
cell2 <- ReadMtx(
  "../analyses/single-cell/input/Pal2021/Exp_data_UMIcounts2.mtx", 
  cells = "../analyses/single-cell/input/Pal2021/Cells2.csv", 
  features = "../analyses/single-cell/input/Pal2021/Genes.txt", 
  cell.column = 1, 
  feature.column = 1,
  cell.sep = ",",
  feature.sep = " ", 
  skip.cell = 1)

# cell 3
cell3 <- ReadMtx(
  "../analyses/single-cell/input/Pal2021/Exp_data_UMIcounts3.mtx", 
  cells = "../analyses/single-cell/input/Pal2021/Cells3.csv", 
  features = "../analyses/single-cell/input/Pal2021/Genes.txt", 
  cell.column = 1, 
  feature.column = 1,
  cell.sep = ",",
  feature.sep = " ", 
  skip.cell = 1)

# cell 4
cell4 <- ReadMtx(
  "../analyses/single-cell/input/Pal2021/Exp_data_UMIcounts4.mtx", 
  cells = "../analyses/single-cell/input/Pal2021/Cells4.csv", 
  features = "../analyses/single-cell/input/Pal2021/Genes.txt", 
  cell.column = 1, 
  feature.column = 1,
  cell.sep = ",",
  feature.sep = " ", 
  skip.cell = 1)

# cell 5
cell5 <- ReadMtx(
  "../analyses/single-cell/input/Pal2021/Exp_data_UMIcounts5.mtx", 
  cells = "../analyses/single-cell/input/Pal2021/Cells5.csv", 
  features = "../analyses/single-cell/input/Pal2021/Genes.txt", 
  cell.column = 1, 
  feature.column = 1,
  cell.sep = ",",
  feature.sep = " ", 
  skip.cell = 1)


# annotation ----

# import annotation tables of different samples in dataset
# samples are organized based on different type of breast cancer

# annotation1
ann1 <- read.csv(
  "../analyses/single-cell/input/Pal2021/Cells1.csv")

# annotation2
ann2 <- read.csv(
  "../analyses/single-cell/input/Pal2021/Cells2.csv")

# annotation3
ann3 <- read.csv(
  "../analyses/single-cell/input/Pal2021/Cells3.csv")

# annotation 4
ann4 <- read.csv(
  "../analyses/single-cell/input/Pal2021/Cells4.csv")

# annotation 5
ann5 <- read.csv(
  "../analyses/single-cell/input/Pal2021/Cells5.csv")


# wrangling ----

# combine annotation
ann_all <- do.call(rbind, list(ann1, ann2, ann3, ann4, ann5))

# type of cells
table(ann_all$cell_type)

# discard cells with no annotation
ann_all_sub <- ann_all[ann_all$cell_type != "", ]

# combine all cells
cell_all <- do.call(cbind, list(cell1, cell2, cell3, cell4, cell5))

# remove single subsets to free up space
rm(cell1, cell2, cell3, cell4, cell5)

# subset all cells keeping only cells with annotation
cell_sub <- cell_all[, colnames(cell_all) %in% ann_all_sub$cell_name]

# check if cell names in ann_all_sub coincide with cell_sub colnames
all(ann_all_sub$cell_name == cell_sub@Dimnames[[2]])


# cell types ids ----

# possible cell types
cell_type <- names(table(ann_all_sub$cell_type))


# split expression ----

# remove cells where all genes == 0
samples_to_keep <- colSums(cell_sub, na.rm = T) != 0

# patients to keep
patient_to_keep <- ann_all_sub$sample[samples_to_keep]
table(patient_to_keep)

# remove genes where all samples == 0
genes_to_keep <- rowSums(cell_sub, na.rm = T) != 0

# susbet exp
exp_sub <- cell_sub[genes_to_keep, samples_to_keep]

rm(cell_sub)
gc()

# transpose
exp_t <- as.data.frame(t(exp_sub))
rm(exp_sub)
gc()

# transform UMI counts into counts per million (CPM)
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

names(prop_nan) <- colnames(exp_pseudo)[-1]

genes_to_keep <- names(prop_nan)[prop_nan < 0.5]

# subset
exp_pseudo <- exp_pseudo[, c("patient", genes_to_keep)]

# save list
saveRDS(
  object = exp_pseudo, 
  file = "/home/ieo5059/mount/local/projects/survival/output/analyses/kaplan-meier/single-cell/output/brca_pal/exp_pseudobulk_all-cells.rds")