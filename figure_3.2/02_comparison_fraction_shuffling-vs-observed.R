# compare fraction of extreme prognosis in low/medium/high expression range
# upon randomization against observed fraction when no shuffling is done
# do it for each cohort


# import packages ----
library(dplyr)
library(ggplot2)


# input variables ----

# survival analysis
surv_analysis <- "DSS"

# cancer cohorts
cohorts <- readRDS("../data/cohorts.rds")

# split by
split_by <- "tertiles"

# significance threshold for multiple testing correction
qvalue_threshold <- 0.2

# method for multiple test correction
method <- "BH"

# significance threshold for log-rank test pairwise comparisons
km_pairwise_threshold <- 0.05

# survival analysis
surv_analysis <- "DSS"

# number of iterations
n_iterations <- 100

# directory with count frequency upon randomization
rand_count_dir <- file.path(
  "../analyses/randomization/count",
  paste0(
    as.character(n_iterations), 
    "_iterations"))


# data ----

# file with frequency of low, medium and high expression range
# with extreme behavior in all cohorts
obs_count <- readRDS(
  "../analyses/counts/new/tertiles/0.2-qvalue-threshold/0.05-km-pairwise-threshold/count_all-cohorts.rds")


# wrangle ----

# exclude cohort where there are no genes associated with extreme prognosis
# in any expression range
cohorts_to_keep <- rowSums(obs_count[, -c(1, 2)]) != 0
obs_count_extreme <- obs_count[cohorts_to_keep, ] 


# loop ----

# loop across cohorts showing at least a gene with extreme prognosis in 
# either low, medium or high expression range
prop_cohort_l <- lapply(seq(nrow(obs_count_extreme)), function(j) {
  
  # cohort
  cohort <- obs_count_extreme$cohort[j]
  
  cat("Process", cohort, "cohort\n") 
  
  # file with N randomization
  rand_count <- readRDS(
    file.path(
      rand_count_dir, 
      cohort, 
      paste0(
        cohort, 
        "_", 
        as.character(n_iterations), 
        "-iterations.rds")))
  
  # subset obs_count to cohort
  obs_count_cohort <- obs_count_extreme[obs_count_extreme$cohort == cohort, ]
  
  # loop over expression range
  freq <- lapply(c("low", "medium", "high"), function(er) {
    
    cat("Look at", er, "expression range\n")
    
    # frequency of observed expression range
    freq_obs_er <- (obs_count_cohort[[paste0("n_", er, "_bad")]] + 
                      obs_count_cohort[[paste0("n_", er, "_good")]])/
      obs_count_cohort$n_genes_tested 
    
    # frequency of background model for expression range
    freq_back_er <- (as.numeric(rand_count[[paste0("n_", er, "_bad")]]) + 
                       as.numeric(rand_count[[paste0("n_", er, "_good")]]))/
      obs_count_cohort$n_genes_tested 
    
    # proportion of times frequency is higher in background model
    res <- mean(freq_back_er > freq_obs_er)
    
    # compute distance between observed frequency and frequencies in background model
    # and normalize by observed value
    d <- abs((freq_obs_er-freq_back_er))/freq_obs_er
    
    # pick the mean of the distance distribution
    # if it is close to 1 it means that for most randomization
    # there are not genes associated with survival by chance
    # the background frequency is 0 and if you subtract 0 to the observed
    # frequency and divide by obs frequency you get 1.
    return(c(res, mean(d)))
    
  })
  
  # frequency as vector
  freq_v <- setNames(
    object = unlist(freq), 
    nm = c(
      "low_proportion_randomization_higher_observed", 
      "low_mean_normalized_distance_observed_background", 
      "medium_proportion_randomization_higher_observed", 
      "medium_mean_normalized_distance_observed_background", 
      "high_proportion_randomization_higher_observed", 
      "high_mean_normalized_distance_observed_background"
    ))
  
  
  # proportion of cases where background frequency > observed frequency
  prop <- c(
    "cohort" = cohort, 
    freq_v)
  
  return(prop)
  
})


# bind elements of list
prop_cohort_df <- as.data.frame(do.call(rbind, prop_cohort_l))

# save
write.table(
  prop_cohort_df, 
  file = file.path(
    "../analyses/randomization/ground_vs_background/frequency",
    as.character(n_iterations), 
    paste0(
      "frequency_ground_vs_background_",
      as.character(n_iterations), 
      ".tsv"
    )
  ), 
  sep = "\t", 
  row.names = F
)


# plot----
ggplot(prop_cohort_df, aes(
  x = cohort, 
  y = medium_proportion_randomization_higher_observed)) + 
  geom_point(stat = "identity", size = 4) + 
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  theme_bw()
