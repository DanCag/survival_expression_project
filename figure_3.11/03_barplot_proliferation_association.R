# barplot
# y axis: number of cohorts where either low, medium or high extreme group
# shows association between survival change and proliferation change 
# between low vs non-low or medium vs non-medium or high vs non-high
# x axis: association (positive/negative/non-significant)


# packages ----
library(ggplot2)

# association low extreme
association_le <- read.delim(
  "../analyses/association/proliferation/new/association_cohort-wise_surv-prol_low.tsv")

# association middle extreme
association_me <- read.delim(
  "../analyses/association/proliferation/new/association_cohort-wise_surv-prol_medium.tsv")

# association high extreme
association_he <- read.delim(
  "../analyses/association/proliferation/new/association_cohort-wise_surv-prol_high.tsv")


# middle extreme ----

# positive association
pa_me <- mean(association_me$pvalue < 0.05 & association_me$odds.ratio > 1)

# negative association
na_me <- mean(association_me$pvalue < 0.05 & association_me$odds.ratio < 1)

# non significant association
ns_me <- mean(association_me$pvalue >= 0.05)


# low extreme ----

# positive association
pa_le <- mean(association_le$pvalue < 0.05 & association_le$odds.ratio > 1)

# negative association
na_le <- mean(association_le$pvalue < 0.05 & association_le$odds.ratio < 1)

# non significant association
ns_le<- mean(association_le$pvalue >= 0.05)


# high extreme ----

# positive association
pa_he <- mean(association_he$pvalue < 0.05 & association_he$odds.ratio > 1)

# negative association
na_he <- mean(association_he$pvalue < 0.05 & association_he$odds.ratio < 1)

# non significant association
ns_he <- mean(association_he$pvalue >= 0.05)


# assemble ----
df <- data.frame(
  "group" = factor(
    c(rep("medium_extreme", 3), rep("low_extreme", 3), rep("high_extreme", 3)),
    levels = c("medium_extreme", "low_extreme", "high_extreme")),
  "association" = factor(
    rep(c("positive", "negative", "non-significant"), 3), 
    levels = c("positive", "negative", "non-significant")),
  "freq_cohort" = c(pa_me, na_me, ns_me, pa_le, na_le, ns_le, pa_he, na_he, ns_he))


# barplot ----
ggplot(
  df, aes(x = association, y = freq_cohort, fill = group)) +
  geom_bar(stat = "identity") +
  facet_wrap(vars(group)) +
  theme_bw()
