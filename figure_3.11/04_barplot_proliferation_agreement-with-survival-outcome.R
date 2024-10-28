# Does the Ki67 expression agrees with the effect of expression range in geneX
# on survival outcome? For instance if "low" expression of geneX is associated
# with better prognosis, does the "non-low" group of patients show an higher 
# expression of Ki67?

# packages ----
library(dplyr)
library(ggplot2)


# parameters ----

# range of interest (either "low", "medium" or "high")
range_interest <- "low"


# input ----

# agreement survival-proliferation
agreement_path <- "../analyses/proliferation/new/shape-wise/0.2-qvalue-threshold/0.05-km-pairwise-threshold/agreement_survival-Ki67-proliferation_genes_expression-range-associated-with-extreme-prognosis.rds"

# association path
association_path <- file.path(
  "../analyses/association/proliferation/new", 
  paste0(
    "association_cohort-wise_surv-prol_", 
    range_interest, 
    ".tsv"))


# data ----

# agreement
agreement <- readRDS(agreement_path)

# association table
association <- read.delim(association_path)


# wrangle ----

# consider only cohorts where the association is significant and positive
cohorts_to_keep <- association$cohort[association$pvalue < 0.05 & association$odds.ratio > 1]

# subset agreement
agreement_sub <- agreement[agreement$cohort %in% cohorts_to_keep, ]


# barplot (range of interest) ----

# keep only those genes where range of interest is involved
agreement_range <- agreement_sub[agreement_sub$range == range_interest, ]

# group by cohort and compute the number and the frequency of genes
# with agreement, disagreement, non significance
agreement_summary <- agreement_range %>%
  group_by(cohort, agreement_surv_prol) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

# transform agreement_surv_prol into a factor
agreement_summary$agreement_surv_prol <- factor(
  agreement_summary$agreement_surv_prol, 
  levels = c(
    "disagree", 
    "Ki67_proliferation_not_significant", 
    "agree")
)

# table with info on number of genes per cohort
info <- agreement_summary %>% 
  group_by(cohort) %>%
  summarise(n = sum(n))

# barplot where x and y axis are flipped
ggplot(data = agreement_summary) +
  geom_bar(
    mapping = aes(x = forcats::fct_rev(cohort),
                  y = freq,
                  fill = agreement_surv_prol),
    stat = "identity",
    position = "stack") +
  scale_fill_manual(values = c("#e63946", "#f1faee", "#457b9d")) +
  labs(
    title = "Frequency agreement survival outcome - Ki67 proliferation",
    subtitle = range_interest,
    y = "frequency",
    x = "cohort") +
  theme_bw() +
  coord_flip() +
  geom_text(
    data = info,
    aes(x = cohort, y = 0, label = n))
