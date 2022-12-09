# Author: Samvida S. Venkatesh
# Date: 30/11/20

library(tidyverse)
library(dplyr)
theme_set(theme_bw())
library(RColorBrewer)

BMI <- read.table("/well/lindgren/UKBIOBANK/samvida/BMI/bmi_primary_care_errors_removed.txt",
                  sep = "\t", header = T)

# Method 1 - logFC from min to max ----
## Plot logFC ----

indivs <- BMI %>% group_by(eid, sex, mean_UKBB_BMI) %>%
  summarise(n_obs = n(),
            min_BMI = min(primarycare_BMI),
            max_BMI = max(primarycare_BMI))

# Subset only those individuals for whom we have multiple measurements
indivs <- subset(indivs, indivs$n_obs > 1)

# Calculate logFC
indivs$logFC <- log2((indivs$max_BMI - indivs$min_BMI) / indivs$min_BMI)

# Group by number of measurements
breakpoints <- c(-Inf, 1, 5, 10, 15, Inf)
names <- c("1", "2-5", "6-10", "11-15", "16+")
indivs$nmeasures_bin <- cut(indivs$n_obs, breaks = breakpoints, labels = names)

# Group by weight class
breakpoints <- c(-Inf, 18.5, 25, 30, 35, 40, Inf)
names <- c("Underweight (< 18.5)", "Normal weight [18.5 - 25)", 
           "Pre-obesity [25 - 30)", "Obesity Class I [30 - 35)", 
           "Obesity Class II [35 - 40)", "Obesity Class III (>= 40)")
indivs$BMI_class <- cut(indivs$mean_UKBB_BMI, breaks = breakpoints, 
                    include.lowest = T, right = F,
                    labels = names)

# Plot distribution of logFC
popn_99 <- quantile(indivs$logFC, 0.99)

pdf("/well/lindgren/UKBIOBANK/samvida/BMI/plots/indiv_outliers/logFC_distribution.pdf")

# WOMEN
p <- subset(indivs, indivs$sex == "F")
# Faceted by number of measures
ggplot(p, aes(x = logFC)) +
  facet_wrap(~nmeasures_bin, scales = "free_y") +
  geom_histogram(position = "identity", fill = "#f8766d") +
  geom_vline(xintercept = popn_99, linetype = "dashed", col = "black") +
  labs(x = "logFC (min to max)", y = "Number of women")
# Faceted by weight class
ggplot(subset(p, !is.na(p$BMI_class)), aes(x = logFC)) +
  facet_wrap(~BMI_class, scales = "free_y") +
  geom_histogram(position = "identity", fill = "#f8766d") +
  geom_vline(xintercept = popn_99, linetype = "dashed", col = "black") +
  labs(x = "logFC (min to max)", y = "Number of women")

# MEN
p <- subset(indivs, indivs$sex == "M")
# Faceted by number of measures
ggplot(p, aes(x = logFC)) +
  facet_wrap(~nmeasures_bin, scales = "free_y") +
  geom_histogram(position = "identity", fill = "#00bfc4") +
  geom_vline(xintercept = popn_99, linetype = "dashed", col = "black") +
  labs(x = "logFC (min to max)", y = "Number of men")
# Faceted by weight class
ggplot(subset(p, !is.na(p$BMI_class)), aes(x = logFC)) +
  facet_wrap(~BMI_class, scales = "free_y") +
  geom_histogram(position = "identity", fill = "#00bfc4") +
  geom_vline(xintercept = popn_99, linetype = "dashed", col = "black") +
  labs(x = "logFC (min to max)", y = "Number of men")

dev.off()

## Calculate thresholds ----

thresholds <- indivs %>% group_by(sex, nmeasures_bin, BMI_class) %>%
  summarise(q99 = quantile(logFC, 0.99)[[1]])
write.table(thresholds, "/well/lindgren/UKBIOBANK/samvida/BMI/indiv_outliers_logFC.txt",
            sep = "\t", row.names = F, quote = F)

# Identify outlier individuals

indivs <- indivs %>% group_by(sex, nmeasures_bin, BMI_class) %>%
  mutate(q99 = quantile(logFC, 0.99)[[1]],
         outlier = logFC > q99)

outliers <- indivs[which(indivs$outlier), "eid"]

## Plot trajectories for outlier individuals ----

BMI <- BMI %>% group_by(eid) %>% mutate(n_obs = n())
# Group by number of measurements
breakpoints <- c(-Inf, 1, 5, 10, 15, Inf)
names <- c("1", "2-5", "6-10", "11-15", "16+")
BMI$nmeasures_bin <- cut(BMI$n_obs, breaks = breakpoints, labels = names)

# Group by weight class
breakpoints <- c(-Inf, 18.5, 25, 30, 35, 40, Inf)
names <- c("Underweight (< 18.5)", "Normal weight [18.5 - 25)", 
           "Pre-obesity [25 - 30)", "Obesity Class I [30 - 35)", 
           "Obesity Class II [35 - 40)", "Obesity Class III (>= 40)")
BMI$BMI_class <- cut(BMI$mean_UKBB_BMI, breaks = breakpoints, 
                          include.lowest = T, right = F,
                          labels = names)

outliers <- subset(BMI, BMI$eid %in% outliers$eid)
outliers$outlier <- T

# Randomly sample from the remaining (non-outlier) population to display trajectories
sample_BMI <- subset(indivs, !indivs$outlier) %>% 
  group_by(sex, nmeasures_bin, BMI_class) %>%
  sample_n(size = 15)
sample_BMI <- subset(BMI, BMI$eid %in% sample_BMI$eid)
sample_BMI$outlier <- F

# Plot trajectories faceted by sex, number of measurements, and weight class

pdata <- bind_rows(outliers, sample_BMI)

plist <- pdata %>% group_by(nmeasures_bin, BMI_class) %>% 
  group_split()

p <- lapply(plist, function (df) {
  bc <- unique(df$BMI_class)
  nb <- unique(df$nmeasures_bin)
  res <- ggplot(subset(df, !df$outlier), 
                aes(x = age_years, y = primarycare_BMI, group = eid,
                    col = sex)) +
    facet_wrap(~sex, nrow = 2, scales = "free") +
    geom_point(col = "#ebebeb") +
    geom_line(col = "#ebebeb", size = 0.7) +
    geom_point(data = subset(df, df$outlier), 
               aes(x = age_years, y = primarycare_BMI, group = eid,
                   col = sex)) +
    geom_line(data = subset(df, df$outlier), 
               aes(x = age_years, y = primarycare_BMI, group = eid,
                   col = sex)) +
    labs(x = "Age (years)", y = "BMI", title = paste("UKBB", bc, 
                                                     "Primary care", nb, "Measurements",
                                                     sep = " "))
  return (res)
})

pdf("/well/lindgren/UKBIOBANK/samvida/BMI/plots/indiv_outliers/seqlogFC_outliers.pdf",
    onefile = T)
print(p)
dev.off()

# Method 2 - Largest sequential logFC ----

# First arrange the measurements by age (they should already be sorted anyway)
# lag from dplyr gets the value in the previous row
seq_change <- BMI %>% group_by(eid) %>%
  arrange(age_years, .by_group = T) %>%
  mutate(seqlogFC = log2(abs(primarycare_BMI - lag(primarycare_BMI)) / 
                           lag(primarycare_BMI)))

# Get max value of change for each individual
indivs <- seq_change %>% group_by(eid, sex, mean_UKBB_BMI) %>% 
  summarise(n_obs = n(),
            max_logFC = max(seqlogFC, na.rm = T))

# Subset only those individuals for whom we were able to calculate logFC
indivs <- subset(indivs, is.finite(indivs$max_logFC))

# Same plots as for logFC (method 1) and all following steps are the same

# Method 3 - >= 2 large sequential logFC ----

# First arrange the measurements by age (they should already be sorted anyway)
# lag from dplyr gets the value in the previous row
seq_change <- BMI %>% group_by(eid) %>%
  arrange(age_years, .by_group = T) %>%
  mutate(seqlogFC = log2(abs(primarycare_BMI - lag(primarycare_BMI)) / 
                           lag(primarycare_BMI)))

# Get every logFC recorded
all_jumps <- seq_change$seqlogFC[is.finite(seq_change$seqlogFC)]

# Summary:
summary(all_jumps)
p95_popn <- quantile(all_jumps, 0.95)[[1]]

# Get number of jumps larger than p95 for each individual
indivs <- seq_change %>% group_by(eid, sex, mean_UKBB_BMI) %>% 
  summarise(n_obs = n(), 
            n_large_jumps = length(which(seqlogFC > p95_popn)))

# Plots same as for methods 1&2, but plot on population level

# Calculate thresholds ----

# Only keep finite logFC
seq_change <- subset(seq_change, is.finite(seq_change$seqlogFC))
# Calculate extreme values
thresholds <- seq_change %>% group_by(sex, nmeasures_bin, BMI_class) %>%
  summarise(q99 = quantile(seqlogFC, 0.99)[[1]])
write.table(thresholds, "/well/lindgren/UKBIOBANK/samvida/BMI/indiv_outliers_seqlogFC.txt",
            sep = "\t", row.names = F, quote = F)

# Identify outlier individuals

indivs <- seq_change %>% group_by(sex, nmeasures_bin, BMI_class) %>%
  mutate(q99 = quantile(seqlogFC, 0.99)[[1]]) 

indivs <- indivs %>% group_by(eid, sex, nmeasures_bin, BMI_class, n_obs) %>%
  summarise(n_large_jumps = length(which(seqlogFC > q99)))

outliers <- indivs[which(indivs$n_large_jumps > 1), "eid"]