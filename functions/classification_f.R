# loop over gene and classify gene expression distributions into: 
# - bimodal
# - normal
# - positively/negatively skewed

# classification function
classification_f <-  function(gene, exp, b_threshold = 0.05, n_threshold = 0.05) {
  
  # cat("Assess", g, "gene distribution\n")
  
  # Bimodality ----
  
  # Hartigans' dip test for unimodality / multimodality
  dip_test <- diptest::dip.test(exp[[gene]])
  
  if (dip_test$p.value < b_threshold) {
    
    # check if distribution has two modes
    is_bimodal <- LaplacesDemon::is.bimodal(exp[[gene]])
    
    if ( is_bimodal ) {
      
      label <- "bimodal"
      pval <- dip_test$p.value
      
    } else {
      
      label <- "unlabeled"
      pval <- NA
      
    }
    
  } else {
    
    # Normality ----
    
    # check for normality
    normality_test <- shapiro.test(exp[[gene]])
    
    if ( normality_test$p.value > n_threshold ) {
      
      # cat("Distribution is normal\n\n")
      
      label <- "normal"
      pval <- normality_test$p.value
      
    } else {
      
      # compute skewness of distribution
      skewness <- moments::skewness(exp[[gene]], na.rm = T)
      
      if (skewness > 1){
        
        # cat("Positive skewness\n\n")
        
        label <- "positive_skewed"
        pval <- NA
        
      } else if (skewness < -1) {
        
        # cat("Negative skewness\n\n")
        
        label <- "negative_skewed"
        pval <- NA
        
      } else if (skewness == 0 ){ # there should not be any gene in this group
        
        # cat("Symmetric\n\n")
        label <- "symmetric"
        pval <- NA
        
      } else {
        
        # cat("Unlabeled\n\n")
        label <- "unlabeled"
        pval <- NA
        
      }
    }
  }
  
  return(c(
    "gene" = gene,
    "label" = label, 
    "p-value" = pval))
  
}
