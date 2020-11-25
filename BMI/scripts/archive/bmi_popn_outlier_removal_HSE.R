# Author: Samvida S. Venkatesh
# Date: 24/11/2020

# Remove extreme BMI values based on HSE (Health Survey England) data
# Values +/- 10% more extreme than the min and max HSE are removed as noise

library(tidyverse)
library(dplyr)
library(RColorBrewer)
theme_set(theme_bw())

# Read data ----

hse_data <- read.table("/well/lindgren/UKBIOBANK/samvida/BMI/HSE_summaries.txt",
                       sep = "\t", header = T)
pc_data <- read.table("/well/lindgren/UKBIOBANK/samvida/BMI/bmi_primary_care_annotated.txt",
                      sep = "\t", header = T)

# Calculate cut-off thresholds ----

# Values +/- 10% more extreme than the min and max HSE 

MIN <- 0.9*min(hse_data$min)
MAX <- 1.1*max(hse_data$max)

# Plot to visualise cut-off thresholds
pdf("/well/lindgren/UKBIOBANK/samvida/BMI/plots/popn_outlier_thresholds.pdf")
ggplot(pc_data, aes(x = primarycare_BMI, color = sex, fill = sex)) +
  geom_density(alpha = 0.25, na.rm = T) +
  geom_vline(xintercept = MIN) +
  geom_vline(xintercept = MAX) +
  geom_vline(xintercept = min(hse_data$min), linetype = "dashed") + 
  geom_vline(xintercept = max(hse_data$max), linetype = "dashed") + 
  labs(x = "BMI", y = "Density") +
  xlim(c(10, 100))

# Remove population-level outliers based on thresholds ----

pc_data$remove <- pc_data$primarycare_BMI < MIN | 
  pc_data$primarycare_BMI > MAX

# What do we lose? ----

removed <- pc_data[pc_data$remove, ]
cleaned <- pc_data[!pc_data$remove, ]

table(removed$sex)

# Look at distribution of BMI observations that were removed for being too low
ggplot(subset(removed, removed$primarycare_BMI < MIN), 
       aes(x = primarycare_BMI, color = sex, fill = sex)) +
  geom_density(alpha = 0.25, na.rm = T) +
  labs(x = "BMI", y = "Density", 
       title = "Distribution of BMI measures removed for being too low") 

# or too high
ggplot(subset(removed, removed$primarycare_BMI > MAX), 
       aes(x = primarycare_BMI, color = sex, fill = sex)) +
  geom_density(alpha = 0.25, na.rm = T) +
  labs(x = "BMI", y = "Density", 
       title = "Distribution of BMI measures removed for being too high") 

cleaned_indivs <- cleaned %>% group_by(eid, sex, mean_UKBB_BMI) %>% 
  summarise(n_obs = n())
removed_indivs_low <- removed[removed$primarycare_BMI <= MIN, ] %>% 
  group_by(eid, sex, mean_UKBB_BMI) %>%
  summarise(n_removed = n())
removed_indivs_high <- removed[removed$primarycare_BMI >= MAX, ] %>% 
  group_by(eid, sex, mean_UKBB_BMI) %>%
  summarise(n_removed = n())

# Look at the UKBIOBANK BMI distribution in individuals who have observation errors
ggplot(cleaned_indivs, aes(x = mean_UKBB_BMI)) +
  geom_density(col = "#1B9E77", fill = "#1B9E77",
               alpha = 0.25, na.rm = T) +
  geom_density(data = removed_indivs_low, 
               col = "#D95F02", fill = "#D95F02",
               alpha = 0.25, na.rm = T) +
  geom_density(data = removed_indivs_high, 
               col = "#7570B3", fill = "#7570B3",
               alpha = 0.25, na.rm = T) +
  labs(x = "BMI", y = "Density") +
  xlim(c(12, 80))
