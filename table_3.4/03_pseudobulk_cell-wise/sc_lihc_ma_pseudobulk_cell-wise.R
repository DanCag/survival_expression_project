# single-cell lihc dataset from Ma et al., 2021 (GSE151530)
# split scRNA-seq into cell contribution
# run this script in the workstation


# packages ----
library(Seurat)
library(dplyr)


# cell ----

# expression sparse matrix
exp <- ReadMtx(
  mtx = "../analyses/single-cell/input/GSE151530/GSE151530_matrix.mtx",
  cells = "../analyses/single-cell/input/GSE151530/GSE151530_barcodes.tsv", 
  features = "../analyses/single-cell/input/GSE151530/GSE151530_genes.tsv", 
  cell.column = 1, 
  feature.column = 2,
  feature.sep = "\t")

# annotation
ann <- read.delim(
  "../analyses/single-cell/input/GSE151530/GSE151530_Info.txt")


# wrangling ----

# keep only cells from hepatocellular carcinoma
ann_hepa <- ann[grep("C", ann$Sample, invert = T), ]

# check the types of cells available
table(ann_hepa$Type)

# remove unclassified cells
ann_hepa_classified <- ann_hepa[ann_hepa$Type != "unclassified", ]

# subset cells in expression
exp_sub <- exp[, colnames(exp) %in% ann_hepa_classified$Cell]

# possible cell types
cell_type <- names(table(ann_hepa_classified$Type))

# list with cell types ids
cell_ids_l <- lapply(cell_type, function(x) {
  
  result <- ann_hepa_classified$Cell[ann_hepa_classified$Type == x]
  
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
  exp_cell <- exp_sub[, colnames(exp_sub) %in% cell_ids_l[[i]]]
  
  # annotation of that particular cell type
  ann_cell <- ann_hepa_classified[ann_hepa_classified$Type == cell_type[i], ]
  
  # check if order of cell names in ann_cell and exp_cell coincides
  all(ann_cell$Cell == colnames(exp_cell))
  
  # transpose
  exp_cell_t <- as.data.frame(t(exp_cell))
  
  # transform UMI counts into counts per million
  exp_cpm <- as.data.frame(t(apply(exp_cell_t, 1, function(x)(x/sum(x))*1e06)))
  
  # add patient column
  exp_cpm <- exp_cpm %>% 
    mutate(
      "patient" = ann_cell$Sample, 
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
  file = "../analyses/single-cell/output/lihc_ma/exp_cell_types_mean-no-zero_l.rds")
