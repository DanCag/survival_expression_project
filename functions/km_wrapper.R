# packages ----
library(survival)
library(survminer)


# Run KM analysis and compute log-rank pvalue
km_pvalue_f <- function(gene, gene_labels, surv_df) {
  
  # add gene labels to surv df
  surv_df[[paste(gene, "cat", sep = "_")]] <- gene_labels
  
  # prepare formula for log-rank test
  f <- as.formula(
    paste(
      "Surv(years, event_int)",
      "~",
      paste(
        gene,
        "cat",
        sep = "_")
    ))
  
  # pairwise comparison
  
  # log-rank test comparison of the three groups
  # pairwise_survdiff() comes from package survminer
  sd <- pairwise_survdiff(
    formula = f,
    data = surv_df)
  
  # transform matrix into a named vector 
  sd_v <- setNames(
    object = as.vector(t(sd$p.value)),
    nm = paste(
      rep(rownames(sd$p.value),
          each = length(colnames(sd$p.value))),
      colnames(sd$p.value),
      sep = "_"))
  
  # discard "low_low", medium_medium", and "high_high" groups
  sd_v <- sd_v[!names(sd_v) %in% c(
    "low_low", "medium_medium", "high_high")]
  
  # fit object  
  fo <- survdiff(
    formula = f, 
    data = surv_df
  )
  
  # extract pvalue
  pv <- fo$pvalue
  
  return(c(pv, sd_v))
  
}


# Run KM analysis and return a fit dataframe object
km_fit_f <- function(gene, gene_labels, surv_df) {
  
  # add gene labels to surv df
  surv_df[[paste(gene, "cat", sep = "_")]] <- gene_labels
  
  
  # prepare formula for log-rank test
  f <- as.formula(
    paste(
      "Surv(years, event_int)",
      "~",
      paste(
        gene,
        "cat",
        sep = "_")
    ))
  
  # define the strata names
  strata_names <- c("low", "medium", "high")
  
  # survfit object
  fit <- survfit(
    formula = f,
    data = surv_df)
  
  # summary object of survfit
  res <- summary(fit)
  
  # desired columns
  cols <- lapply(c(2:6, 8:11) , function(x) res[x])
  
  # survift object into a df
  fit_df <- do.call(data.frame, cols)
  
  # strata order
  strata_order <- paste0(
    gene,
    "=",
    strata_names)
  
  # order fit_df based on strata column levels
  fit_df_ordered <- fit_df %>%
    arrange(factor(strata, levels = strata_order))
  
  return(fit_df_ordered)
  
}