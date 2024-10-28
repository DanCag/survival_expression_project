# single-cell brain dataset from Neftel et al. 2019
# split scRNA-seq into cell contribution
# run this script in the workstation


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
  "../analyses/single-cell/input/Neftel2019/Cells.csv"
)


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

# list with cell types ids
cell_ids_l <- lapply(cell_type, function(x) {
  
  result <- ann$cell_name[ann$cell_type == x]
  
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
  exp_cell <- exp[, colnames(exp) %in% cell_ids_l[[i]]]
  
  # transpose
  exp_cell_t <- as.data.frame(t(exp_cell))
  
  # remove cells where all genes == 0
  samples_to_keep <- rowSums(exp_cell_t, na.rm = T) != 0
  
  # subset exp_cell_t
  exp_cell_sub_t <- exp_cell_t[samples_to_keep, ]
  
  # annotation of that particular cell type
  ann_cell <- ann[ann$cell_type %in% cell_type[i], ]
  
  # subset annotation keeping only samples_to_keep == True
  ann_cell_sub <- ann_cell[ann_cell$cell_name %in% names(samples_to_keep)[samples_to_keep], ]
  
  # check if order of cell names in ann_cell and exp_cell coincides
  all(ann_cell_sub$cell_name == rownames(exp_cell_sub_t))
  
  # add patient column
  exp_cell <- exp_cell_sub_t %>% 
    mutate(
      "patient" = ann_cell_sub$sample, 
      .before = colnames(exp_cell_sub_t)[1])

  # # group by patient and normalize considering only non 0 entries
  exp_pseudo <- exp_cell %>% 
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
  
})

# assign names to elements of list
exp_df_cell_type_l <- setNames(
  object = exp_df_cell_type_l, 
  nm = cell_type
)

# save list
saveRDS(
  object = exp_df_cell_type_l, 
  file = "../analyses/single-cell/output/gbm_neftel/exp_cell_types_mean-no-zero_l.rds")
