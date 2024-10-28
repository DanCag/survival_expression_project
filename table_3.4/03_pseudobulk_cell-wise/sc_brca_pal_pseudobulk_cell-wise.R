# single-cell breast dataset from Pal et al. 2021
# split scRNA-seq into cell contribution

# run this script in workstation


# packages ----
library(Seurat)


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

# list with cell types ids
cell_ids_l <- lapply(cell_type, function(x) {
  
  result <- ann_all_sub$cell_name[ann_all_sub$cell_type == x]
  
})

# assign names to cell_ids_l
cell_ids_l <- setNames(
  object = cell_ids_l, 
  nm = cell_type)


# split expression ----

# loop over cell types and split expression matrix
exp_df_cell_type_l <- lapply(seq(cell_ids_l), function(i) {
  
  cat("Process", cell_type[i], "\n")
  
  # expression of that particular cell type
  exp_sub <- cell_sub[, colnames(cell_sub) %in% cell_ids_l[[i]]]
  
  # transpose
  exp_sub_t <- as.data.frame(t(exp_sub))
  
  # transform in CPM
  cpm <- as.data.frame(t(apply(exp_sub_t, 1, function(x)((x/sum(x)) * 1e06))))
  
  # add patient column
  cpm <- cpm %>% 
    mutate(
      "patient" = sub("_[^_]+$", "", rownames(cpm)), 
      .before = colnames(cpm)[1])
  
  # exp_sub_t[1:4, 1:4]
  
  cpm_pseudo <- cpm %>% 
    group_by(patient) %>% 
    summarise_all(~mean(.[. > 0]))
  
  # proportion of NaN
  prop_nan <- sapply(
    seq(ncol(cpm_pseudo))[-1], 
    function(j) {mean(is.nan(cpm_pseudo[[j]]))})
  
  names(prop_nan) <- colnames(cpm_pseudo)[-1]
  
  # # keep genes with less than 50 % dropouts
  # genes_to_keep2 <- names(prop_zero)[prop_zero < 0.5]
  genes_to_keep <- names(prop_nan)[prop_nan < 0.5]
  
  # subset exp_patient
  cpm_pseudo_sub <- cpm_pseudo[, c("patient", genes_to_keep)]

  
})

# assign names to elements of list
exp_df_cell_type_l <- setNames(
  object = exp_df_cell_type_l, 
  nm = cell_type
)

# save list
saveRDS(
  object = exp_df_cell_type_l, 
  file = "../analyses/single-cell/output/brca_pal/exp_cell_types_mean-no-zero_l.rds")