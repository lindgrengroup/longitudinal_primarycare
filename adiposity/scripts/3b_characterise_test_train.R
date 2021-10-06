# Author: Samvida S. Venkatesh
# Date: 05/10/21

library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

# Read test set, training set, and covariates ----

test_adipo <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/TEST_SET_adiposity.rds")
train_adipo <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/TRAINING_SET_adiposity.rds")
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/covariates.rds")

PHENOTYPES <- names(covars)

# Build test vs training characteristics tables -----

demo_table <- lapply(PHENOTYPES, function (p) {
  test_df <- subset(covars[[p]], covars[[p]]$eid %in% test_adipo[[p]]$eid)
  test_df$cut <- "test"
  
  train_df <- subset(covars[[p]], covars[[p]]$eid %in% train_adipo[[p]]$eid)
  train_df$cut <- "training"
  
  df <- bind_rows(test_df, train_df)
  
  # Get mean and S.D. of numeric columns
  summary_df <- df %>% group_by(cut) %>% 
    summarise(across(c(baseline_age, baseline_trait, FUyrs, FU_n), 
                     list(mean = mean, sd = sd))) 
  
  # Calculate percentage female
  perc_female <- df %>% group_by(cut, sex) %>%
    summarise(n = n()) %>%
    mutate(perc = n/sum(n)) %>% 
    pivot_wider(names_from = sex,
                values_from = c(n, perc))
  summary_df <- merge(summary_df, perc_female, by = "cut")
  
  # From full adiposity data
  test_full <- test_adipo[[p]]
  test_full$cut <- "test"
  train_full <- train_adipo[[p]]
  train_full$cut <- "training"
  all <- bind_rows(test_full, train_full)
  # Calculate percentage of each data provider
  perc_dp <- all %>% group_by(cut, data_provider) %>%
    summarise(n = n()) %>%
    mutate(perc = n/sum(n)) %>% 
    pivot_wider(names_from = data_provider,
                values_from = c(n, perc))
  summary_df <- merge(summary_df, perc_dp, by = "cut")
  
  write.table(summary_df, 
              paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/test_train_characteristics_", 
                     p, ".txt"),
              sep = "\t", quote = F, row.names = F)
  return (summary_df)
})

# Plot summaries of test vs training characteristics ----

## Sex split ----

sex_splits <- lapply(PHENOTYPES, function (p) {
  test_df <- subset(covars[[p]], covars[[p]]$eid %in% test_adipo[[p]]$eid)
  test_df$cut <- "test"
  train_df <- subset(covars[[p]], covars[[p]]$eid %in% train_adipo[[p]]$eid)
  train_df$cut <- "training"
  df <- bind_rows(test_df, train_df)
  
  # Calculate percentage by sex
  perc_female <- df %>% group_by(cut, sex) %>%
    summarise(n = n()) %>%
    mutate(perc = n/sum(n)) 
  
  perc_female$biomarker <- p
  return (perc_female)
})
sex_splits <- bind_rows(sex_splits)

sex_plot <- ggplot(sex_splits, aes(x = sex, y = perc,
                                    fill = cut, color = cut)) +
  facet_wrap(~biomarker, scales = "free") +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Sex", y = "Fraction", title = "Sex split in test vs training")

## Data provider split ----

dp_splits <- lapply(PHENOTYPES, function (p) {
  # From full adiposity data
  test_full <- test_adipo[[p]]
  test_full$cut <- "test"
  train_full <- train_adipo[[p]]
  train_full$cut <- "training"
  all <- bind_rows(test_full, train_full)
  # Calculate percentage of each data provider
  
  perc_dp <- all %>% group_by(cut, data_provider) %>%
    summarise(n = n()) %>%
    mutate(perc = n/sum(n)) 
  
  perc_dp$biomarker <- p
  return (perc_dp)
})
dp_splits <- bind_rows(dp_splits)

dp_plot <- ggplot(dp_splits, aes(x = data_provider, y = perc,
                                   fill = cut, color = cut)) +
  facet_wrap(~biomarker, scales = "free") +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Data provider", y = "Fraction", 
       title = "Data provider split in test vs training")

## Numerical baseline variables ----

covar_comparisons_df <- lapply(PHENOTYPES, function (p) {
  # Create joined df for plots
  test_df <- subset(covars[[p]], covars[[p]]$eid %in% test_adipo[[p]]$eid)
  test_df$cut <- "test"
  train_df <- subset(covars[[p]], covars[[p]]$eid %in% train_adipo[[p]]$eid)
  train_df$cut <- "training"
  df <- bind_rows(test_df, train_df)
  df$biomarker <- p
  return (df)
})
covar_comparisons_df <- bind_rows(covar_comparisons_df)

# Baseline age
bl_age_plot <- ggplot(covar_comparisons_df, 
                      aes(x = cut, y = baseline_age)) +
  facet_wrap(~biomarker, scales = "free") +
  geom_violin(aes(fill = cut), position = position_dodge(1)) +
  geom_boxplot(width = 0.1) + 
  scale_y_continuous() +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Cut", y = "Baseline age") +
  theme(legend.position = "none")

# Baseline trait value
bl_trait_plot <- ggplot(covar_comparisons_df, 
                      aes(x = cut, y = baseline_trait)) +
  facet_wrap(~biomarker, scales = "free") +
  geom_violin(aes(fill = cut), position = position_dodge(1)) +
  geom_boxplot(width = 0.1) + 
  scale_y_continuous() +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Cut", y = "Baseline trait") +
  theme(legend.position = "none")

# Number of FU years
fuyr_plot <- ggplot(covar_comparisons_df, aes(x = cut, y = FUyrs, fill = cut)) +
  facet_wrap(~biomarker, scales = "free") +
  geom_boxplot(position = position_dodge(1)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Cut", y = "Follow-up length (years)") +
  theme(legend.position = "none")

# Number of FU measurements
n_fu_plot <- ggplot(covar_comparisons_df, aes(x = cut, y = FU_n, fill = cut)) +
  facet_wrap(~biomarker, scales = "free") +
  geom_boxplot(position = position_dodge(1)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Cut", y = "# follow-up measures") +
  theme(legend.position = "none")

# Print all plots to a pdf ----

pdf("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/test_train_characteristics.pdf",
    onefile = T)
print(sex_plot)
print(dp_plot)
print(bl_age_plot)
print(bl_trait_plot)
print(fuyr_plot)
print(n_fu_plot)
dev.off()
