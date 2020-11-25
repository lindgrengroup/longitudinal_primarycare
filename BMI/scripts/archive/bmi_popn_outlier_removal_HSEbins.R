## ARCHIVED ##
# Author: Samvida S. Venkatesh
# Date: 18/11/2020

# Remove extreme BMI values based on HSE (Health Survey England) data
# Values +/- 10% more extreme than the min and max HSE in each age group and sex
# are removed as measurement error / noise

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

# Values +/- 10% more extreme than the min and max HSE in each 
# age group and sex (across the 1993 and 2018 surveys)

thresholds <- hse_data %>% group_by(sex, age_bin) %>% 
  summarise(min_cutoff = 0.9*min(min),
            max_cutoff = 1.1*max(max))

# Sort primary care data into the same bins as HSE data ----

breakpoints <- c(16, 25, 35, 45, 55, 65, 75, Inf)
names <- c("16-24", "25-34", "35-44", "45-54", "55-64",
           "65-74", "75+")
pc_data$age_bin <- cut(pc_data$age_years, 
                       include.lowest = T, right = F,
                       breaks = breakpoints, labels = names)

# Remove population-level outliers based on thresholds ----

remove_idx <- c()
for (i in 1:dim(thresholds)[1]) {
  idx <- which(pc_data$sex == thresholds$sex[i] & 
                 pc_data$age_bin == thresholds$age_bin[i])
  tmp <- pc_data[idx, ]
  rmv <- which(tmp$primarycare_BMI < thresholds$min_cutoff[i] | 
                 tmp$primarycare_BMI > thresholds$max_cutoff[i])
  remove_idx <- c(remove_idx, idx[rmv])
}

pc_data$remove <- F
pc_data$remove[remove_idx] <- T

# Summary statistics of data to be kept vs removed ----

summ <- pc_data %>% group_by(sex, age_bin, remove) %>%
  summarise(n = n(), 
            min_BMI = min(primarycare_BMI), max_BMI = max(primarycare_BMI),
            mean_BMI = mean(primarycare_BMI), 
            median_BMI = median(primarycare_BMI))

write.table(summ, 
            "/well/lindgren/UKBIOBANK/samvida/BMI/plots/popn_outlier_summary.txt", 
            sep = "\t", quote = F, row.names = F)

# Number of outlier measurements per individual ----

p <- pc_data %>% group_by(eid, sex, mean_UKBB_BMI) %>%
  summarise(n_obs = n(),
            n_outliers = sum(remove == T)) %>%
  mutate(ratio_outliers_obs = n_outliers / n_obs)

# Obviously some individuals are getting excluded that should not be
# CHANGE THRESHOLDING 