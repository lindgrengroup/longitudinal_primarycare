# Author: Samvida S. Venkatesh
# Date: 09/12/20

library(tidyverse)
library(dplyr)
theme_set(theme_bw())
library(RColorBrewer)

# Read and supplement data ----

BMI <- read.table("/well/lindgren/UKBIOBANK/samvida/BMI/bmi_primary_care_errors_removed.txt",
                  sep = "\t", header = T)

# Number of BMI measurements
BMI <- BMI %>% group_by(eid) %>% mutate(n_obs = n())

# Group by number of measurements
breakpoints <- c(-Inf, 1, 5, 10, 15, Inf)
names <- c("1", "2-5", "6-10", "11-15", "16+")
BMI$nobs_bin <- cut(BMI$n_obs, breaks = breakpoints, labels = names)

# Group by mean UKBIOBANK BMI obesity class
breakpoints <- c(-Inf, 18.5, 25, 30, 35, 40, Inf)
names <- c("Underweight (< 18.5)", "Normal weight [18.5 - 25)", 
           "Pre-obesity [25 - 30)", "Obesity Class I [30 - 35)", 
           "Obesity Class II [35 - 40)", "Obesity Class III (>= 40)")
BMI$BMI_class <- cut(BMI$mean_UKBB_BMI, breaks = breakpoints, 
                     include.lowest = T, right = F,
                     labels = names)

# Calculate sequential logFC
seq_change <- BMI %>% group_by(eid) %>%
  arrange(age_years, .by_group = T) %>%
  mutate(FC = abs(primarycare_BMI - lag(primarycare_BMI)) / 
           lag(primarycare_BMI),
         age_change = age_years - lag(age_years),
         logFC = log2(FC),
         seqlogFC = log2(FC / age_change))

# Calculate outliers ----

# Only keep values where the FC is not NA (NA: first obs in the series or when
# only one observation available)
seq_change <- subset(seq_change, !is.na(seq_change$FC))

# 99th percentile of individuals in each category:

summary(seq_change$logFC)
popn_99_logFC <- quantile(seq_change$logFC, 0.99)[[1]]

summary(seq_change$seqlogFC)
# Don't look at positive infinite values (as these come from age-changes of 0)
finseqlogFC <- seq_change$seqlogFC[which(seq_change$seqlogFC != Inf)]
popn_99_seqlogFC <- quantile(finseqlogFC, 0.99, na.rm = T)[[1]]

seq_change$outlier <- seq_change$logFC > popn_99_logFC &
  seq_change$seqlogFC > popn_99_seqlogFC

outliers <- unique(seq_change$eid[which(seq_change$outlier)])
outliers <- subset(BMI, BMI$eid %in% outliers)
outliers$outlier <- T

# Plot trajectories for outliers -----

# Randomly sample from the remaining (non-outlier) population in each
# sex, nobs_bin, and obesity class to display trajectories
sample_BMI <- subset(seq_change, !seq_change$outlier) %>% 
  distinct(eid, sex, nobs_bin, BMI_class) %>%
  group_by(sex, nobs_bin, BMI_class) %>%
  sample_n(size = min(n(), 50))
sample_BMI <- subset(BMI, BMI$eid %in% sample_BMI$eid)
sample_BMI$outlier <- F

# Plot trajectories faceted by sex and weight class

pdata <- bind_rows(outliers, sample_BMI)
plist <- pdata %>% group_by(nobs_bin, BMI_class) %>% group_split()

# For lines showing UKBIOBANK obesity class thresholds
mins <- c(-Inf, 18.5, 25, 30, 35, 40)
maxes <- c(18.5, 25, 30, 35, 40, Inf)
names <- c("Underweight (< 18.5)", "Normal weight [18.5 - 25)", 
           "Pre-obesity [25 - 30)", "Obesity Class I [30 - 35)", 
           "Obesity Class II [35 - 40)", "Obesity Class III (>= 40)")
names(mins) <- names
names(maxes) <- names

p <- lapply(plist, function (df) {
  bc <- unique(df$BMI_class)
  nb <- unique(df$nobs_bin)
  res <- ggplot(subset(df, !df$outlier), 
                aes(x = age_years, y = primarycare_BMI, group = eid,
                    col = sex)) +
    facet_wrap(~sex, nrow = 2, scales = "free") +
    geom_point(col = "#ebebeb") +
    geom_line(col = "#ebebeb", size = 0.7) +
    geom_hline(yintercept = mins[bc], linetype = "dashed") +
    geom_hline(yintercept = maxes[bc], linetype = "dashed") +
    geom_point(data = subset(df, df$outlier), 
               aes(x = age_years, y = primarycare_BMI, group = eid,
                   col = sex)) +
    geom_line(data = subset(df, df$outlier), 
              aes(x = age_years, y = primarycare_BMI, group = eid,
                  col = sex)) +
    labs(x = "Age (years)", y = "BMI", title = paste("UKBB", bc, 
                                                     "Primary care", nb, 
                                                     "Measurements",
                                                     sep = " "))
  return (res)
})

pdf("/well/lindgren/UKBIOBANK/samvida/BMI/plots/indiv_outliers/combined_outliers.pdf",
    onefile = T)
print(p)
dev.off()

# Remove noisy observation for individuals identified as outliers ----

outliers <- outliers %>% group_by(eid) %>% 
  mutate(med_primarycare = median(primarycare_BMI),
         dist_to_med = abs(primarycare_BMI - med_primarycare)) 
# If there were only two observations, calculate the distance to mean UKBB BMI
outliers$dist_to_med[which(outliers$n_obs < 3)] <- 
  abs(outliers$primarycare_BMI[which(outliers$n_obs < 3)] - 
        outliers$mean_UKBB_BMI[which(outliers$n_obs < 3)])

outliers <- outliers %>% mutate(remove = dist_to_med == max(dist_to_med))

cleaned <- merge(BMI, outliers[, c("eid", "event_dt", "primarycare_BMI",
                                   "remove")], 
                 all.x = T)
cleaned$remove[is.na(cleaned$remove)] <- F
cleaned <- subset(cleaned, !cleaned$remove)

# Repeat the above steps with the cleaned BMI data for a second round ----

BMI2 <- cleaned[, c(1:7)]

# Number of BMI measurements
BMI2 <- BMI2 %>% group_by(eid) %>% mutate(n_obs = n())

# Group by number of measurements
breakpoints <- c(-Inf, 1, 5, 10, 15, Inf)
names <- c("1", "2-5", "6-10", "11-15", "16+")
BMI2$nobs_bin <- cut(BMI2$n_obs, breaks = breakpoints, labels = names)

# Group by mean UKBIOBANK BMI obesity class
breakpoints <- c(-Inf, 18.5, 25, 30, 35, 40, Inf)
names <- c("Underweight (< 18.5)", "Normal weight [18.5 - 25)", 
           "Pre-obesity [25 - 30)", "Obesity Class I [30 - 35)", 
           "Obesity Class II [35 - 40)", "Obesity Class III (>= 40)")
BMI2$BMI_class <- cut(BMI2$mean_UKBB_BMI, breaks = breakpoints, 
                     include.lowest = T, right = F,
                     labels = names)

# Plot trajectories in each sex, obesity class, number of measures 
# in order to find the last few noisy values:

plist <- BMI2 %>% group_by(nobs_bin, BMI_class) %>% group_split()

# For lines showing UKBIOBANK obesity class thresholds
mins <- c(-Inf, 18.5, 25, 30, 35, 40)
maxes <- c(18.5, 25, 30, 35, 40, Inf)
names <- c("Underweight (< 18.5)", "Normal weight [18.5 - 25)", 
           "Pre-obesity [25 - 30)", "Obesity Class I [30 - 35)", 
           "Obesity Class II [35 - 40)", "Obesity Class III (>= 40)")
names(mins) <- names
names(maxes) <- names

p <- lapply(plist, function (df) {
  bc <- unique(df$BMI_class)
  nb <- unique(df$nobs_bin)
  res <- ggplot(df, aes(x = age_years, y = primarycare_BMI, group = eid,
                    col = sex)) +
    facet_wrap(~sex, nrow = 2, scales = "free") +
    geom_point(col = "#ebebeb") +
    geom_line(col = "#ebebeb", size = 0.7) +
    geom_hline(yintercept = mins[bc], linetype = "dashed") +
    geom_hline(yintercept = maxes[bc], linetype = "dashed") +
    labs(x = "Age (years)", y = "BMI", title = paste("UKBB", bc, 
                                                     "Primary care", nb, 
                                                     "Measurements",
                                                     sep = " "))
  return (res)
})

pdf("/well/lindgren/UKBIOBANK/samvida/BMI/plots/indiv_outliers/trajectories.pdf",
    onefile = T)
print(p)
dev.off()


