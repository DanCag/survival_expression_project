# single-cell breast dataset from Wu et al. 2021
# split scRNA-seq into cell contribution
# run this script in the workstation


# packages ----
library(Seurat)
library(dplyr)


# cell ----

# expression sparse matrix
exp <- ReadMtx(
  mtx = "../analyses/single-cell/input/Wu2021/Exp_data_UMIcounts.mtx",
  cells = "../analyses/single-cell/input/Wu2021/Cells.csv", 
  features = "../analyses/single-cell/input/Wu2021/Genes.txt", 
  cell.column = 1, 
  cell.sep = ",",
  skip.cell = 1,
  feature.column = 1)

# metadata
meta <- read.csv(
  "../analyses/single-cell/input/Wu2021/Meta-data.csv")

# annotation
ann <- read.csv(
  "../analyses/single-cell/input/Wu2021/Cells.csv"
)


# metadata ----

# cancer type
table(meta$cancer_type)

# number of patients
length(table(meta$patient))


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
  
  # annotation of that particular cell type
  ann_cell <- ann[ann$cell_type == cell_type[i], ]
  
  # check if order of cell names in ann_cell and exp_cell coincides
  all(ann_cell$Cell == colnames(exp_cell))
  
  # transpose
  exp_cell_t <- as.data.frame(t(exp_cell))
  
  # transform UMI counts into counts per million
  exp_cpm <- as.data.frame(t(apply(exp_cell_t, 1, function(x)(x/sum(x))*1e06)))
  
  # add patient column
  exp_cpm <- exp_cpm %>% 
    mutate(
      "patient" = ann_cell$patient, 
      .before = colnames(exp_cpm)[1])
  
  exp_pseudo <- exp_cpm %>% 
    group_by(patient) %>% 
    summarise_all(~mean(.[. > 0]))
  
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
  file = "../analyses/single-cell/output/brca_wu/exp_cell_types_mean-no-zero_l.rds")
