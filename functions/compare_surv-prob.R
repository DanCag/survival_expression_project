# ++++++++++++++++++++ #
# Survival probability #
# ++++++++++++++++++++ #

# packages ----
library(dplyr)

# compare survival probability of fit object
# for different stratification at a certain time point (year = 5 by default)
# go backwards in time if there is not info at year 5
# a difference <= 0.5 year is tolerated when comparing survival of different strata
# strata has 3 groups (low, medium, high)
surv_prob_comp_3_groups <- function(fit_df, cg, year = 5, strata_names) {
  
  # look for suitable time to compare survival probabilities
  # of high, medium and low group
  while (year > 1) { # if around year there are not suitable
                     # times (high, medium and low) to compare, 
                     # year = year - 1
                     # I set the minimum year to look at = 2 because at year = 1
                     # curves might not separate yet
    
    # compute times
    
    # first group
    time_first <- fit_df %>% 
      group_by(strata) %>% 
      slice(which.min(abs(time - year))) %>% 
      filter(strata == paste0(cg, "=", strata_names[1])) %>% 
      pull(time)
    
    # second group
    time_second <- fit_df %>% 
      group_by(strata) %>% 
      slice(which.min(abs(time - year))) %>% 
      filter(strata == paste0(cg, "=", strata_names[2])) %>% 
      pull(time) 
    
    # third group
    time_third <- fit_df %>% 
      group_by(strata) %>% 
      slice(which.min(abs(time - year))) %>% 
      filter(strata == paste0(cg, "=", strata_names[3])) %>% 
      pull(time) 
    
    # if there is not suitable time (empty objects)
    if ( length(time_first) == 0 |
         length(time_second) == 0 |
         length(time_third) == 0 ){
      
      # decrease year by one
      year <- year - 1
      
      # assign NAs to times and survival probabilities
      time_first <- NA
      time_second <- NA
      time_third <- NA
      surv_first <- NA
      surv_second <- NA
      surv_third <- NA
      
    } else {
      
      # check if times are close enough one with the other
      if ( abs(time_first - time_second) <= 0.5 &
           abs(time_first - time_third) <= 0.5 &
           abs(time_second - time_third) <= 0.5 ) {
        
        # compute survivals
        
        # first
        surv_first <- fit_df %>% 
          group_by(strata) %>% 
          slice(which.min(abs(time - year))) %>% 
          filter(strata == paste0(cg, "=", strata_names[1])) %>% 
          pull(surv)
        
        # second
        surv_second <- fit_df %>% 
          group_by(strata) %>% 
          slice(which.min(abs(time - year))) %>% 
          filter(strata == paste0(cg, "=", strata_names[2])) %>% 
          pull(surv) 
        
        # third
        surv_third <- fit_df %>% 
          group_by(strata) %>% 
          slice(which.min(abs(time - year))) %>% 
          filter(strata == paste0(cg, "=", strata_names[3])) %>% 
          pull(surv) 
        
        # if found suitable times and survival, 
        # stop while loop
        break
        
      } else {
        
        # decrease by 1 year
        year <- year - 1
        
        # assign NAs to times and survival probabilities
        time_first <- NA
        time_second <- NA
        time_third <- NA
        surv_first <- NA
        surv_second <- NA
        surv_third <- NA
      }
    }
    
  }
  
  # times
  times <- setNames(
    object = c(time_first, time_second, time_third), 
    nm = paste0("time_", strata_names)
  )
  
  # survival probability
  sp <- setNames(
    object = c(surv_first, surv_second, surv_third), 
    nm = paste0("surv_", strata_names)
  )
  
  # list time - survival probability
  time_sp_l <- setNames(
    object = list(times, sp), 
    nm = c("time", "sp"))
  
  return(time_sp_l)
}