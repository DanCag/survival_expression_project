# GO enrichment analysis on cohorts 
# genes where medium expression range is associated with extreme prognosis

# packages ----
library(clusterProfiler)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(dplyr)


# parameters ----

# survival analysis
surv_analysis <- "DSS"

# are you interested both at bad/good prognosis
# or at one of the two only?
# value can be "good", "bad" or "both"
prog_interest <- "both"


# input ----

# directory with pvalue_info
pvalue_dir <- "../analyses/p-value/new/shape-wise/all-genes"


# data ----

# dataframe with significant genes obtained upon shape-wise split
sig_genes <- readRDS(
  "../analyses/genes_list/new/shape-wise/significant/extreme/significant-genes_shape-wise.rds")


# wrangle ----

# keep only medium extreme genes
sig_genes_medium <- sig_genes[
  grep(paste(c("low", "high"), collapse = "|"),
       sig_genes$prognosis, invert = T), ]

if ( prog_interest == "bad" ) {
  
  # keep only bad prognosis genes
  sig_genes_sub <- sig_genes_medium[grep("bad", sig_genes_medium$prognosis), ]
  
} else if ( prog_interest == "good" ) {
  
  # keep only good prognosis genes
  sig_genes_sub <- sig_genes_medium[grep("good", sig_genes_medium$prognosis), ]
  
} else {
  
  # keep both bad and good prognosis
  sig_genes_sub <- sig_genes_medium
  
}

# number of middle extreme (bad/good/both) genes per cohort
table(sig_genes_sub$cohort)

# consider only cohorts with at leas one middle extreme gene
cohorts <- names(table(sig_genes_medium$cohort))[
  table(sig_genes_medium$cohort) >= 1]


# list with GO enrichment for each cohort
ego_middle_l <- lapply(seq(cohorts), function(i) {
  
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
      ".rds")
  )
  
  # check if file exists for that specific cohort
  if ( file.exists(pvalue_path) ) {
    
    # import datasets ----
    
    # import pvalue of all genes tested
    pvalue <- readRDS(pvalue_path)
    
    # genes tested
    genes_tested <- pvalue$gene
    
    # keep middle extreme genes of that cohort
    middle_extreme_genes <- sig_genes_sub$gene[
      sig_genes_sub$cohort == cohorts[i]]
    
    # GO enrichment
    ego <- enrichGO(
      gene = middle_extreme_genes, 
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
names(ego_middle_l) <- cohorts

# bind
ego_middle_df <- bind_rows(ego_middle_l)

# save enrichment result
saveRDS(
  ego_middle_df,
  file = "../analyses/enrichment/cluster_profile/enriched-terms_medium-extreme.rds")

# number of occurence of each term
table(ego_middle_df$ID)
# no shared terms

# save list of GO ids to feed REVIGO
write.table(
  ego_middle_df$ID, 
  file = "../analyses/enrichment/cluster_profile/go_ids_enriched_medium-extreme.txt", 
  quote = F, 
  row.names = F, 
  col.names = F)
