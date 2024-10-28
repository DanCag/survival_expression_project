# single-cell hnsc dataset from Puram et al. 2020
# split scRNA-seq into cell contribution
# run this script in the workstation

# work on it 

# Process expression of HNSC patients 

# packages ----
library(dplyr)


# data ----

# expression + annotation
all <- read.delim(
  gzfile("../analyses/single-cell/input/GSE103322_HNSCC_all_data.txt.gz"))


# annotation ----

# separate annotation from gene expression
ann <- as.data.frame(t(all[1:5, -1]))

# assign colnames
colnames(ann) <- gsub(" ", "_", all[1:5, 1])

ann[1:5, ]
table(ann$`non-cancer_cell_type`)

# correct -Fibroblast into Fibroblast
ann$`non-cancer_cell_type` <- ifelse(
  ann$`non-cancer_cell_type` == "-Fibroblast", "Fibroblast", 
  ann$`non-cancer_cell_type`)
table(ann$`non-cancer_cell_type`)

# transform cancer cell type column
ann$`non-cancer_cell_type` <- ifelse(
  ann$`non-cancer_cell_type` == 0, "cancer", ann$`non-cancer_cell_type`)


# expression ----

# separate expression
exp <- all[6:nrow(all), ]
colnames(exp)[1] <- "SYMBOL"
exp$SYMBOL <- gsub("'", "", exp$SYMBOL)
exp[1:4, 1:5]

# convert column into numeric
exp[, 2:ncol(exp)] <- apply(
  exp[, 2:ncol(exp)], 2, function(x) as.numeric(as.character(x)))

str(exp)


# wrangling ----

# cells from lymph nodes
lymph_node_cells <- rownames(ann)[ann$Lymph_node == 1]

# remove lymph node cells from exp
exp_no_lymph <- exp[, !colnames(exp) %in% lymph_node_cells]
exp_no_lymph[1:4, 1:4]

# cell types
cell_type <- names(table(ann$`non-cancer_cell_type`))

# list with cell types ids
cell_ids_l <- lapply(cell_type, function(x) {
  
  result <- rownames(ann)[ann$`non-cancer_cell_type` == x]
  
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
  exp_sub <- exp_no_lymph[, colnames(exp_no_lymph) %in% c("SYMBOL", cell_ids_l[[i]])]
  
  # transpose
  exp_sub_t <- as.data.frame(t(exp_sub[, -1]))
  colnames(exp_sub_t) <- exp_sub$SYMBOL
  exp_sub_t[1:4, 1:4]
  
  # add patient column
  exp_sub_t <- exp_sub_t %>% 
    mutate("patient" = sub("_.*", "", rownames(exp_sub_t)), .before = colnames(exp_sub_t)[1])
  
  # transform HNSCC into HN
  exp_sub_t$patient <- sub("SCC", "", exp_sub_t$patient)
  
  exp_pseudo <- exp_sub_t %>% 
    group_by(patient) %>% 
    summarise_all(~mean(.[. > 0]))
  
  # proportion of NaN
  prop_nan <- sapply(
    seq(ncol(exp_pseudo))[-1], 
    function(j) {mean(is.nan(exp_pseudo[[j]]))})
  
  names(prop_nan) <- colnames(exp_pseudo)[-1]
  
  # # keep genes with less than 50 % dropouts
  genes_to_keep <- names(prop_nan)[prop_nan < 0.5]
  
  # subset exp_patient
  exp_pseudo_sub <- exp_pseudo[, c("patient", genes_to_keep)]
  
})

# assign names to elements of list
exp_df_cell_type_l <- setNames(
  object = exp_df_cell_type_l, 
  nm = cell_type
)

# save list
saveRDS(
  object = exp_df_cell_type_l, 
  file = "../analyses/single-cell/output/hnsc_puram/exp_cell_types_mean-no-zero_l.rds")
