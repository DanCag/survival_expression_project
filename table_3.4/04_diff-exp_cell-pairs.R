# compare gene expression distributions between cell types

# import packages ----
library(dplyr)


# parameter ----

# cancer type
cancer_type <- "skcm_li"


# data ----

# list with expression dataframe of different cell types
exp_df_cell_type_l <- readRDS(
  file.path(
    "../analyses/single-cell/output", 
    cancer_type, 
    "exp_cell_types_mean-no-zero_l.rds")
)


# wrangling ----

# genes expressed in cell types
genes_l <- lapply(names(exp_df_cell_type_l), function(ct) {
  
  cat("Process", ct, "\n")
  
  # dataframe
  df <- exp_df_cell_type_l[[ct]]
  
  # genes 
  genes <- colnames(df)[-1]
  
})

# genes with less than 50% 0s
genes_to_keep <- Reduce(intersect, genes_l)


# differential expression new ----

# list of all possible cell types pairs
comb_l <- combn(names(exp_df_cell_type_l), m = 2, simplify = F)

# loop over cell types
diff_exp_df_l <- lapply(seq(comb_l), function(i) {
  
  cat("Compare", comb_l[[i]][1], comb_l[[i]][2], "\n")
  
  # datasets
  df1 <- as.data.frame(exp_df_cell_type_l[[comb_l[[i]][1]]])
  df2 <- as.data.frame(exp_df_cell_type_l[[comb_l[[i]][2]]])
  
  # loop over genes and run wilcoxon
  diff_exp_l <- lapply(seq(genes_to_keep), function(j){
    
    # gene
    gene <- genes_to_keep[j]
    
    cat("Compare", gene, j, "\n")
    
    # wilcoxon test
    wilcoxon <- wilcox.test(df1[[gene]], df2[[gene]])
    
    # if wilcoxon is significant
    if ( wilcoxon$p.value < 0.05 ) {
      
      # compute means
      mean1 <- mean(df1[[gene]], na.rm = T)
      mean2 <- mean(df2[[gene]], na.rm = T)
      
      # compare means
      if ( mean1 > mean2 ) {
        
        res <- "higher"
        
      } else {
        
        res <- "lower"
      }
      
      info <- c(
        gene,
        "yes", 
        res, 
        wilcoxon$p.value)
      
    } else {
      
      info <- c(
        gene, 
        "no", 
        NA, 
        wilcoxon$p.value)
      
    }
    
  })
  
  # dataframe with result of wilcoxon test
  diff_exp_df <- as.data.frame(do.call(rbind, diff_exp_l))
  
  # assign names to columns
  colnames(diff_exp_df) <- c(
    "gene", 
    "significant", 
    paste0(comb_l[[i]][1], "_vs_", comb_l[[i]][2]), 
    "pvalue")
  
  # transform pvalue column into numeric
  diff_exp_df$pvalue <- as.numeric(diff_exp_df$pvalue)
  
  # add padjust column
  diff_exp_df <- diff_exp_df %>% 
    mutate("p_adjust_bh" = p.adjust(method = "BH", diff_exp_df$pvalue), .after = "pvalue")
  
  # add padjust bonferroni
  diff_exp_df <- diff_exp_df %>% 
    mutate("p_adjust_bonf" = p.adjust(method = "bonferroni", diff_exp_df$pvalue), .after = p_adjust_bh)
  
  return(diff_exp_df)
})

# assign names to elements of list
diff_exp_df_l <- setNames(
  object = diff_exp_df_l,
  nm = paste0(sapply(comb_l, "[[", 1), "_", sapply(comb_l, "[[", 2)))

# save
saveRDS(
  object = diff_exp_df_l,
  file = file.path(
    "../output/analyses/kaplan-meier/single-cell/output", 
    cancer_type, 
    "diff_exp_df_mean-no-zero_l.rds"))