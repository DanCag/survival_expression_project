# Prcess single-cell HNSC dataset from Puram et al., 2017 (GSE103322)
# run this script in new workstation (more powerful than old workstation)

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

# remove cells where all genes == 0
samples_to_keep <- colSums(exp_no_lymph[, -1], na.rm = T) != 0

# patients to keep
patient_to_keep <- rownames(ann)[rownames(ann) %in% names(samples_to_keep)]
table(patient_to_keep)

# remove genes where all samples == 0
genes_to_keep <- rowSums(exp_no_lymph[, -1], na.rm = T) != 0

# subset exp
exp_sub <- exp_no_lymph[genes_to_keep, samples_to_keep]


# transpose
exp_t <- as.data.frame(t(exp_sub[, -1]))
colnames(exp_t) <- exp_sub[[1]]
# rm(exp_sub)
# gc()

# add patient column
exp_t <- exp_t %>% 
  mutate("patient" = sub("_.*", "", rownames(exp_t)), .before = colnames(exp_t)[1])

# transform HNSCC into HN
exp_t$patient <- sub("SCC", "", exp_t$patient)

# group by patient and normalize considering only non 0 entries
exp_pseudo <- exp_t %>% 
  group_by(patient) %>% 
  summarise_all(~mean(.[. > 0]))

# proportion of NaN
prop_nan <- sapply(
  seq(ncol(exp_pseudo))[-1], 
  function(j) {mean(is.nan(exp_pseudo[[j]]))})

# assign gene names to prop_nan
names(prop_nan) <- colnames(exp_pseudo)[-1]

# # keep genes with less than 50 % NaN
genes_to_keep <- names(prop_nan)[prop_nan < 0.5]

# subset exp_patient
exp_pseudo_sub <- exp_pseudo[, c("patient", genes_to_keep)]

# save list
saveRDS(
  object = exp_pseudo_sub, 
  file = "../analyses/single-cell/output/hnsc_puram/exp_pseudobulk_all-cells.rds")
