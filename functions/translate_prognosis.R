# return string with prognosis info
# according to which expression range shows extreme prognosis
translate_prognosis <- function(count_df, gene) {
  
  # subset count_df to that gene
  count_sub <- count_df[count_df$gene == gene, ]
  
  # prognosis 
  if ( count_sub$count_low == 1) {
    
    if (count_sub$count_medium == -1) {
      
      prog <- "low_good_medium_bad"
      
    } else if (count_sub$count_high == -1) {
      
      prog <- "low_good_high_bad"
      
    } else {
      
      prog <- "low_good"
    }
    
  } else if ( count_sub$count_low == -1 ) {
    
    if (count_sub$count_medium == 1) {
      
      prog <- "low_bad_medium_good"
      
    } else if (count_sub$count_high == 1) {
      
      prog <- "low_bad_high_good"
      
    } else {
      prog <- "low_bad"
    }
    
  } else { # low is 0
    
    if (count_sub$count_medium == -1) {
      
      if (count_sub$count_high == 1) {
        
        prog <- "medium_bad_high_good"
        
      } else {
        
        prog <- "medium_bad"
        
      }
      
    } else if (count_sub$count_medium == 1) {
      
      if (count_sub$count_high == -1) {
        
        prog <- "medium_good_high_bad"
        
      } else {
        
        prog <- "medium_good"
        
      }
      
    } else if (count_sub$count_high == 1 & count_sub$count_medium == 0) {
      
      prog <- "high_good"
      
    } else if (count_sub$count_high == -1 & count_sub$count_medium == 0) {
      
      prog <- "high_bad"
      
    }
    
  }
  
  return(prog)
  
}

