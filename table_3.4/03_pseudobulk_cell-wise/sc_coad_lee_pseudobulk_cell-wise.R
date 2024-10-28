# single-cell coad dataset from Lee et al, 2020 (GSE132465)
# split expression of CRC patients into expression of different cell types
# run this script in the workstation

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

# remove tumor cells
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

# list with cell types ids
cell_ids_l <- lapply(cell_type, function(x) {
  
  result <- ann$Index[ann$Cell_type == x]
  
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
  exp_cell <- as.data.frame(tpm_t[rownames(tpm_t) %in% cell_ids_l[[i]], ])
  
  # annotation of that particular cell type
  ann_cell <- ann_tumor[ann_tumor$Index %in% cell_ids_l[[i]], ]
  
  # check if exp_cell cells are in the same order as ann_cell
  all(rownames(exp_cell) == ann_cell$Index)
  
  # add patient column
  exp_cell <- exp_cell %>% 
    mutate("patient" = ann_cell$Patient, .before = colnames(exp_cell)[1])
  
  exp_pseudo <- exp_cell %>% 
    group_by(patient) %>% 
    summarise_all(~mean(.[. > 0]))
  
  # proportion of NaN
  prop_nan <- sapply(
    seq(ncol(exp_pseudo))[-1], 
    function(j) {mean(is.nan(exp_pseudo[[j]]))})
  
  names(prop_nan) <- colnames(exp_pseudo)[-1]
  
  # # keep genes with less than 50 % dropouts
  # genes_to_keep2 <- names(prop_zero)[prop_zero < 0.5]
  genes_to_keep <- names(prop_nan)[prop_nan < 0.5]
  
  # subset exp_patient
  exp_pseudo_sub <- exp_pseudo[, c("patient", genes_to_keep)]
  
  
})

# assign names to elements of list
exp_df_cell_type_l <- setNames(
  object = exp_df_cell_type_l, 
  nm = cell_type
)

# remove mast cells because there are too few genes
dim(exp_df_cell_type_l$`Mast cells`)
exp_df_cell_type_l$`Mast cells` <- NULL

# save list
saveRDS(
  object = exp_df_cell_type_l, 
  file = "../analyses/single-cell/output/coad_lee/exp_cell_types_mean-no-zero_l.rds")
