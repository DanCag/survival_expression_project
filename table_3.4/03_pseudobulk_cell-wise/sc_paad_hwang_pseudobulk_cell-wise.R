# single-cell pdca dataset from Hwang et al. 2022
# split scRNA-seq into cell contribution
# run this script in the workstation


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

# list with cell types ids
cell_ids_l <- lapply(cell_type, function(x) {
  
  result <- ann_all$cell_name[ann_all$cell_type == x]
  
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
  exp_sub <- cell_all[, colnames(cell_all) %in% cell_ids_l[[i]]]
  
  # transpose
  exp_sub_t <- as.data.frame(t(exp_sub))
  
  # remove cells where all genes == 0
  samples_to_keep <- rowSums(exp_sub_t, na.rm = T) != 0
  
  all(samples_to_keep)
  
  # ann cell
  ann_cell <- ann_all[ann_all$cell_type == cell_type[i], ]
  
  all(ann_cell$cell_name == rownames(exp_sub_t))
  
  # add patient column
  tp10k <- exp_sub_t %>% 
    mutate(
      "patient" = ann_cell$sample, 
      .before = colnames(exp_sub_t)[1])
  
  # exp_sub_t[1:4, 1:4]
  
  tp10k_pseudo <- tp10k %>% 
    group_by(patient) %>% 
    summarise_all(~mean(.[. > 0]))
  
   # proportion of NaN
  prop_nan <- sapply(
    seq(ncol(tp10k_pseudo))[-1], 
    function(j) {mean(is.nan(tp10k_pseudo[[j]]))})
  
  names(prop_nan) <- colnames(tp10k_pseudo)[-1]
  
  # # keep genes with less than 50 % dropouts
  genes_to_keep <- names(prop_nan)[prop_nan < 0.5]
  
  # subset exp_patient
  tp10k_pseudo_sub <- tp10k_pseudo[, c("patient", genes_to_keep)]
  

})

# assign names to elements of list
exp_df_cell_type_l <- setNames(
  object = exp_df_cell_type_l, 
  nm = cell_type
)

# remove some cell types (< 2000 genes detected)
exp_df_cell_type_l$B_cell <- NULL
exp_df_cell_type_l$Mast <- NULL
exp_df_cell_type_l$Neutrophil <- NULL
exp_df_cell_type_l$Plasma <- NULL

# save list
saveRDS(
  object = exp_df_cell_type_l, 
  file = "../analyses/single-cell/output/paad_hwang/exp_cell_types_mean-no-zero_l.rds")
