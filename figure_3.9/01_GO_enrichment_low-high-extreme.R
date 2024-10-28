# GO enrichment analysis on cohorts 
# genes where either low or high expression range (peripheral) is associated 
# with extreme prognosis


# packages ----
library(clusterProfiler)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(dplyr)
library(ggplot2)


# parameters ----

# survival analysis
surv_analysis <- "DSS"


# input ----

# directory with pvalue_info
pvalue_dir <- "../analyses/p-value/new/shape-wise/all-genes"


# data ----

# dataframe with significant genes obtained upon shape-wise split
sig_genes <- readRDS(
  "../analyses/genes_list/new/shape-wise/significant/extreme/significant-genes_shape-wise.rds")


# wrangle ----

# keep only low/high extreme genes
sig_genes_sub <- sig_genes[
  grep("medium", sig_genes$prognosis, invert = T), ]

# number of genes per cohort
table(sig_genes_sub$cohort)

# consider only cohorts with at least one peripheral extreme gene
cohorts <- names(table(sig_genes_sub$cohort))[
  table(sig_genes_sub$cohort) >= 1]


# list with GO enrichment for each cohorts
ego_l <- lapply(seq(cohorts), function(i) {
  
  cat("Process", cohorts[i], "\n")
  
  # p-value path
  pvalue_path <- file.path(
    pvalue_dir, 
    cohorts[i], 
    surv_analysis, 
    paste0(
      "p-value_", 
      cohorts[i],
      "_", 
      surv_analysis, 
      ".rds"))
  
  # check if file exists for that specific cohort
  if ( file.exists(pvalue_path) ) {
    
    # import datasets ----
    
    # import pvalue of all genes tested
    pvalue <- readRDS(pvalue_path)
    
    # genes tested
    genes_tested <- pvalue$gene
    
    # keep extreme genes of that cohort
    extreme_genes <- sig_genes_sub$gene[
      sig_genes_sub$cohort == cohorts[i]]
    
    # GO enrichment
    ego <- enrichGO(
      gene = extreme_genes, 
      universe = genes_tested, 
      keyType = "ENSEMBL", 
      ont = "ALL",
      OrgDb = organism,
      pvalueCutoff = 0.01)
    
    if ( !is.null(ego) ) {
      
      # extract table with over representation analysis
      ego_table <- ego@result
      
      if ( nrow(ego_table) > 0 ) {
        
        # simplify terms
        ego <- enrichplot::pairwise_termsim(ego)
        ego_simplified <- simplify(
          ego, 
          cutoff = 0.7, 
          by = "p.adjust", 
          select_fun = min)
        
        # table with over representation analysis
        ego_table_simplified <- ego_simplified@result
        
        # add column with cohort info repeated
        ego_table_simplified <- ego_table_simplified %>% 
          mutate(cohort = rep(cohorts[i], nrow(ego_table_simplified)),
                 .before = ONTOLOGY)
        
        # result
        res <- ego_table_simplified
        
      } else {
        
        cat("No enriched terms\n")
        res <- NULL
      }
      
    } else {
      
      cat("GO enrichment is NULL")
      res <- NULL
      
    }
    
  } else {
    
    cat("Files needed are not there\n")
    
    res <- NULL
  }
  
  return(res)
})

# assign cohorts' names to elements of the list
names(ego_l) <- cohorts

# bind
ego_df <- bind_rows(ego_l)

# save ego_df
saveRDS(
  ego_df, 
  file = "../analyses/enrichment/cluster_profile/enriched-terms_low-high-extreme.rds")

# save list of GO ids to feed REVIGO
write.table(
  ego_df$ID, 
  file = "../analyses/enrichment/cluster_profile/go_ids_enriched_low-high-extreme.txt", 
  quote = F, 
  row.names = F, 
  col.names = F)

# number of occurrence of each term
table(ego_df$ID)
shared_go <- names(table(ego_df$ID))[table(ego_df$ID) > 1]

# keep only those terms that are shared by at least two cohorts
ego_shared_df <- ego_df[ego_df$ID %in% shared_go, ]

# extract GO description and number of occurrence
df <- data.frame(
  "GO_id" = as.character(names(table(ego_shared_df$ID))), 
  "id_description" = as.character(names(table(ego_shared_df$Description))), 
  "count" = as.numeric(table(ego_shared_df$Description)))

# order by count
df <- df %>% 
  arrange(desc(count))

# save list of GO ids to feed REVIGO
write.table(
  df$GO_id, 
  file = "../analyses/enrichment/cluster_profile/go_ids_shared_enriched_low-high-extreme.txt", 
  quote = F, 
  row.names = F, 
  col.names = F)
