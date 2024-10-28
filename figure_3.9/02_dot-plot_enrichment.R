# barplot enrichment
# either for peripheral extreme genes (figure 3.9a)
# or middle extreme genes (figure 3.9b)


# packages ----
library(dplyr)
library(ggplot2)


# parameters ----

# what are you looking at? peripheral or middle extrem
exp_range <- "peripheral"


# data ----

# dataframe with enrichment results for all cohorts 

if (exp_range == "peripheral") {

  ego_df <- readRDS(
  "../output/analyses/kaplan-meier/enrichment/cluster_profile/enriched-terms_low-high-extreme.rds")
  
} else {
  
  ego_df <- readRDS(
    "../output/analyses/kaplan-meier/enrichment/cluster_profile/enriched-terms_medium-extreme.rds")
  
}

# plot ----

# number of occurrence of each term
table(ego_df$ID)

# share go

if (exp_range == "peripheral"){
  
  shared_go <- names(table(ego_df$ID))[table(ego_df$ID) > 3]
  
} else {
  
  shared_go <- names(table(ego_df$ID))[table(ego_df$ID) > 0]
  
}

# keep only those terms that are shared by at least two cohorts
ego_shared_df <- ego_df[ego_df$ID %in% shared_go, ]

# remove terms with less than 10 genes
ego_shared_sub_df <- ego_shared_df[ego_shared_df$Count >= 5, ]

# number of genes of GO IDs
ego_shared_sub_df$n_genes <- as.numeric(
  sapply(strsplit(x = ego_shared_sub_df$GeneRatio, split = "/"), "[[", 2))

# extract GO description and average number of occurrence and p.adjusted across
# cohorts where the term pops up
df <- ego_shared_sub_df %>% 
  group_by(ID, Description) %>% 
  summarize(
    count_average = round(mean(Count)),
    n_average = round(mean(n_genes)),
    pvalue_average = mean(p.adjust))

# compute ratio
df$gene_ratio <- df$count_average/df$n_average

# order by pvalue
df_ordered <- df[order(df$pvalue_average, decreasing = F), ]
# df_ordered <- df[order(df$gene_ratio, decreasing = T), ]

# subset only top 15 terms in term of count
df_top15 <- df_ordered[1:15, ]

# remova NAs
df_top15 <- na.omit(df_top15)


# barplot
ggplot(df_top15) +
  geom_point(
    aes(
      x = reorder(Description, +gene_ratio), 
      y = gene_ratio, 
      size = count_average, 
      color = pvalue_average),
    stat = "identity") +
  ylim(0, 0.1) +
  coord_flip() + 
  labs(x = "GO term description") + 
  theme_bw()

# dim for saving 
# width: 525
# height: 500
