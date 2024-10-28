# +++++++++++++++++++ #
# Prognostic behavior #
# +++++++++++++++++++ #

# For each gene we provide in input a list with time
# and survival probability info.

# Based on survival probability values, we assign which expression group
# ("low", "medium", "high") is associated with extreme prognosis (bad/good). 
# In case one expression group has an extreme prognostic behavior, 
# and is significantly different from the other two 
# a +1/-1 is assigned to the gene.
# If the other two groups are not significantly different, 
# then they are not considered to be associated with extreme prognosis.

# tertiles split taking into account distribution shapes
counts_3_groups <- function(time_sp_l, pairwise_info, strata_names, threshold){
  
  if ( all(is.na(time_sp_l[["sp"]])) ) { # - if survivals are NAs  
    
    # assign counts
    counts <- setNames(
      object = rep(NA, 3), 
      nm = paste0("count_", strata_names)
    )
    
  } else {  # - if survivals are not NAs 
    
    # 3 possible scenarios:
    # 1. second group is significantly different from both first and second
    # 2. third group is significantly different from both second and first
    # 3. first group is significantly different from both second and third
    
    # 1. second group is significantly different from both first and third 
    if (pairwise_info[names(pairwise_info) == "medium_low"] < threshold &
        pairwise_info[names(pairwise_info) == "high_medium"] < threshold) {
      
      ## comparisons of survival (first, second, third)
      
      # - if second is the lowest
      if (
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] < 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] & 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] < 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"]
      ) {
        
        # check if third and first comparison is significant          
        if (pairwise_info[names(pairwise_info) == "high_low"] < threshold) {
          
          # comparison of survival (high, low)
          # 3 possible scenarios:
          # a) third > first
          # b) first > third
          # c) third = first (very unlikely) 
          #    this means that even if first and second are significantly different,
          #    at that specific time the two survival probabilities coincides 
          #    and we cannot discriminate which has the more extreme behavior
          
          # a)
          if ( time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] >
               time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] ) { 
            
            # assign counts
            counts <- setNames(
              object = c(0, -1, 1), 
              nm = paste0("count_", strata_names))
            
            # b)
          } else if ( time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] <
                      time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] ) {
            
            # assign counts
            counts <- setNames(
              object = c(1, -1, 0), 
              nm = paste0("count_", strata_names))
            
            # c)  
          } else if ( time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] == 
                      time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] ){
            
            # assign counts
            counts <- setNames(
              object = c("to_disentangle", -1, "to_disentangle"), 
              nm = paste0("count_", strata_names))
            
          }
          
        } else { # third and first are not significantly different
                 # then only second has an extreme behavior
          
          # assign counts
          counts <- setNames(
            object = c(0, -1, 0), 
            nm = paste0("count_", strata_names))
        }
        
      } else if ( # - second is the highest
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] > 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] & 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] >
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] ) {
        
        # check if first and third comparison is significant
        if ( 
          pairwise_info[names(pairwise_info) == "high_low"] < threshold ) {
          
          # comparison of survival (high, low)
          if ( 
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] > 
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] ) { # - if first > third
            
            # assign counts
            counts <- setNames(
              object = c(0, 1, -1), 
              nm = paste0("count_", strata_names))
            
          } else if ( 
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] <
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] ){ # - if first < third
            
            # assign counts
            counts <- setNames(
              object = c(-1, 1, 0), 
              nm = paste0("count_", strata_names))
            
          } else if (
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] == 
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] ){ 
            
            # first == third at that specific time, but later they separate
            
            # assign counts
            counts <- setNames(
              object = c("to_disentangle", 1, "to_disentangle"), 
              nm = paste0("count_", strata_names))
            
          }
          
        } else { # first and third are not significantly different
                 # then only second group has an extreme behavior
          
          # assign counts
          counts <- setNames(
            object = c(0, 1, 0), 
            nm = paste0("count_", strata_names))
          
        }
        
        
      } else if ( # - second is in the middle and third has best prognosis
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] > 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] & 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] < 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] ) { 
        
        # no need to check pvalue between first and third, 
        # because if second group is in the middle
        # and is significantly different from first and third, 
        # then first and third will be also significant
        
        # assign count
        counts <- setNames(
          object = c(-1, 0, 1), 
          nm = paste0("count_", strata_names))
        
      } else if ( # - second is in the middle and first has best prognosis
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] > 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] &
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] < 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] ){ 
        
        # assign count
        counts <- setNames(
          object = c(1, 0, -1), 
          nm = paste0("count_", strata_names))
        
      } else if (# - second is equal to first
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] ==
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"]) {
        
        if ( # second > third
          time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] > 
          time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] ) {
          
          # assign count
          counts <- setNames(
            object = c("to_disentangle", "to_disentangle", -1), 
            nm = paste0("count_", strata_names))
          
        } else if ( # second < third
          time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] < 
          time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] ) {
          
          # assign count
          counts <- setNames(
            object =  c("to_disentangle", "to_disentangle", 1), 
            nm = paste0("count_", strata_names))
          
        }
      } else if ( # second is equal to third
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] == 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] ) {
        
        if ( # if second > first
          time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] > 
          time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"]) {
          
          # assign count
          counts <- setNames(
            object = c(-1, "to_disentangle", "to_disentangle"), 
            nm = paste0("count_", strata_names))
          
        } else if ( # second < first
          time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] < 
          time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"]) {
          
          # assign count
          counts <- setNames(
            object = c(1, "to_disentangle", "to_disentangle"), 
            nm = paste0("count_", strata_names))
        }
        
      }
      
      # 2. third is significantly different from second and first  
    } else if (  
      pairwise_info[names(pairwise_info) == "high_medium"] < threshold &
      pairwise_info[names(pairwise_info) == "high_low"] < threshold ) {
      
      ## comparisons of survival (first, second, third)
      
      # - if third is the lowest
      if (
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] < 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] & 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] <
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"]) {
        
        # check if medium and low comparison is significant          
        if (
          pairwise_info[names(pairwise_info) == "medium_low"] < threshold ) {
          
          # comparison of survival (first, second)
          # 3 possible scenarios:
          # a) first > second
          # b) second > first
          # c) first = second
          
          # a)
          if (
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] > 
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"]) { 
            
            # assign counts
            counts <- setNames(
              object = c(1, 0, -1), 
              nm = paste0("count_", strata_names))
            
            # b)  
          } else if (
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] < 
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"]) { 
            
            # assign counts
            counts <- setNames(
              object = c(0, 1, -1), 
              nm = paste0("count_", strata_names))
            
            # c)  
          } else if (
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] == 
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"]) {
            
            # they are equal at that time, but later the two curves separate
            
            # assign counts
            counts <- setNames(
              object = c("to_disentangle", "to_disentangle", -1), 
              nm = paste0("count_", strata_names))
            
          }
          
        } else { # second and first are not significantly different
                 # then only third has an extreme behavior
          
          # assign counts
          counts <- setNames(
            object = c(0, 0, -1), 
            nm = paste0("count_", strata_names))
        }
        
      } else if ( # - third is the highest
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] > 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] & 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] > 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] ) {
        
        # check if second and first comparison is significant
        if (
          pairwise_info[names(pairwise_info) == "medium_low"] < threshold) {
          
          # comparison of survival (third, first)
          # 2 possible scenarios:
          # a) first > second
          # b) second > first
          # c) second = first
          
          # a)
          if (
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] > 
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] ) {
            
            # assign counts
            counts <- setNames(
              object = c(0, -1, 1), 
              nm = paste0("count_", strata_names))
            
            # b)  
          } else if (
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] <
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] ){ 
            
            # assign counts
            counts <- setNames(
              object = c(-1, 0, 1), 
              nm = paste0("count_", strata_names))
            
            # c) 
          } else if (
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] == 
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"]) {
            # they are equal at that time, but later the two curves separate
            
            # assign counts
            counts <- setNames(
              object = c("to_disentangle", "to_disentangle", 1), 
              nm = paste0("count_", strata_names))
          }
          
        } else { # second and first are not significantly different
                 # then only third has an extreme behavior
          
          # assign counts
          counts <- setNames(
            object = c(0, 0, 1), 
            nm = paste0("count_", strata_names))
          
        }
        
      } else if ( # - third is in the middle and second has best prognosis
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] > 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] & 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] < 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] ) { 
        
        # no need to check pvalue between first and second, 
        # because if third curve is in the middle
        # and is significantly different from first and second, 
        # then first and second will also be significant
        
        # assign count
        counts <- setNames(
          object = c(-1, 1, 0), 
          nm = paste0("count_", strata_names))
        
      } else if ( # - third is in the middle and first has best prognosis
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] > 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] &
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] < 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] ){ 
        
        # assign count
        counts <- setNames(
          object = c(1, -1, 0), 
          nm = paste0("count_", strata_names))
        
      } else if ( # third is identical to first at that specific time 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] ==
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] ) {
        
        if ( # if second > third
          time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] > 
          time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"]) {
          
          # assign count
          counts <- setNames(
            object = c("to_disentangle", 1, "to_disentangle"), 
            nm = paste0("count_", strata_names))
          
        } else if ( # second < third
          time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] < 
          time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"]) {
          
          # assign count
          counts <- setNames(
            object = c("to_disentangle", -1, "to_disentangle"), 
            nm = paste0("count_", strata_names))
        }
        
      } else if ( # third is identical to second at that specific time
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] ==
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"]
      ) {
        
        if ( # if first > third
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] > 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"]) {
        
        # assign count
        counts <- setNames(
          object = c(1, "to_disentangle", "to_disentangle"), 
          nm = paste0("count_", strata_names))
        
      } else if ( # first < third
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] < 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"]) {
        
        # assign count
        counts <- setNames(
          object = c(-1, "to_disentangle", "to_disentangle"), 
          nm = paste0("count_", strata_names)) }
      }
      
      # 3. first is significantly different from third and second  
    } else if (
      pairwise_info[names(pairwise_info) == "high_low"] < threshold &
      pairwise_info[names(pairwise_info) == "medium_low"] < threshold ) {
      
      ## comparisons of survival (first, second, third)
      
      # - if first is the lowest
      if (
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] < 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] & 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] <
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] ) {
        
        # check if second and third comparison is significant          
        if (
          pairwise_info[names(pairwise_info) == "high_medium"] < threshold) {
          
          # comparison of survival (second, third)
          # 3 possible scenarios:
          # a) third > second
          # b) second > third
          # c) second = third (at that specific time, curves separate later)
          
          # a)
          if (
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] > 
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"]) { 
            
            # assign counts
            counts <- setNames(
              object = c(-1, 0, 1), 
              nm = paste0("count_", strata_names))
            
            # b)    
          } else if (
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] <
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] ) {
            
            # assign counts
            counts <- setNames(
              object = c(-1, 1, 0), 
              nm = paste0("count_", strata_names))
            
            # c)
          } else if ( # third is identical to second at that specific time
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] == 
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] ) {
            
            # assign counts
            counts <- setNames(
              object = c(-1, "to_disentangle", "to_disentangle"), 
              nm = paste0("count_", strata_names))
          }
          
        } else { # second and third are not significantly different
                 # then only first has an extreme behavior
          
          # assign counts
          counts <- setNames(
            object = c(-1, 0, 0), 
            nm = paste0("count_", strata_names))
        }
        
      } else if ( # - first is the highest
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] >
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] & 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] > 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] ) {
        
        # check if medium and high comparison is significant
        if (
          pairwise_info[names(pairwise_info) == "high_medium"] < threshold) {
          
          # comparison of survival (second, third)
          # 2 possible scenarios:
          # a) third > second
          # b) second > third
          # c) second = third
          
          # a)
          if (
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] > 
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"]) {
            
            # assign counts
            counts <- setNames(
              object = c(1, -1, 0), 
              nm = paste0("count_", strata_names))
            
            # b)  
          } else if (
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] < 
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] ){
            
            # assign counts
            counts <- setNames(
              object = c(1, 0, -1), 
              nm = paste0("count_", strata_names))
            
            # c)   
          } else if (
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] == 
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] ) {
            
            # assign counts
            counts <- setNames(
              object = c(1, "to_disentangle", "to_disentangle"), 
              nm = paste0("count_", strata_names))
            
          }
          
        } else { # second and third are not significantly different
          # then only first has an extreme behavior
          
          # assign counts
          counts <- setNames(
            object = c(1, 0, 0), 
            nm = paste0("count_", strata_names))
        }
        
        
      } else if ( # - first is in the middle and second has best prognosis
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] > 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] & 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] < 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] ) { 
        
        # no need to check pvalue between second and third, 
        # because if first curve is in the middle
        # and is significantly different from second and third, 
        # then second and third will also be significant
        
        # assign count
        counts <- setNames(
          object = c(0, 1, -1), 
          nm = paste0("count_", strata_names))
        
      } else if ( # - first is in the middle and third has best prognosis
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] > 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] &
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] <
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] ){ 
        
        # assign count
        counts <- setNames(
          object = c(0, -1, 1), 
          nm = paste0("count_", strata_names))
        
      } else if ( # first is significantly different from both second and third
                  # but at that time point its survival is identical to 
                  # second and to disentangle the question we need to plot km 
        (time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] == 
         time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"]) 
      ) {
        
        if ( # if third > first
          time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] > 
          time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] ) {
          
          # assign count
          counts <- setNames(
            object = c("to_disentangle", "to_disentangle", 1), 
            nm = paste0("count_", strata_names))
          
        } else if (# if third < first
          time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] < 
          time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"]) {
          
          # assign count
          counts <- setNames(
            object = c("to_disentangle", "to_disentangle", -1), 
            nm = paste0("count_", strata_names))
        }

      } else if ( # first is significantly different from both second and third
        # but at that time point its survival is identical to 
        # third and to disentangle the question we need to plot km 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] == 
        time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_high"] 
        ) {
          
          if ( # if second > first
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] > 
            time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"] ) {
              
              # assign count
              counts <- setNames(
                object = c("to_disentangle", 1, "to_disentangle"), 
                nm = paste0("count_", strata_names))
              
            } else if (# if second < first
              time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_medium"] < 
              time_sp_l[["sp"]][names(time_sp_l$sp) == "surv_low"]) {
              
              # assign count
              counts <- setNames(
                object = c("to_disentangle", -1, "to_disentangle"), 
                nm = paste0("count_", strata_names))
            }
          
        }
      
    } else { # only one comparison is significant
             # for instance low-high is significant, 
             # but low-medium and medium-high are not significant
      
      # assign count
      counts <- setNames(
        object = c(0, 0, 0), 
        nm = paste0("count_", strata_names))
      
    }
  }
  
  # return vector with times and counts
  return(c(time_sp_l[["time"]],
           time_sp_l[["sp"]],
           counts))
}