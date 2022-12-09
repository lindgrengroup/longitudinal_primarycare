# Author: Samvida S. Venkatesh
# Date: 26/03/21

library(tidyverse)
library(lubridate)
library(lme4)
theme_set(theme_bw())

# Read files ----

adiposity <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/visually_QCd_adiposity.rds")
adiposity <- lapply(adiposity, function (x) {
  x$eid <- as.character(x$eid)
  return (x)
})

general_covars <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/QCd_demographic_covariates.txt",
                             sep = "\t", header = T, comment.char = "$",
                             stringsAsFactors = F)
general_covars <- subset(general_covars, general_covars$sex == "M" & 
                           general_covars$ancestry == "white")

PHENOTYPES <- names(adiposity)
p <- "WHR"

# Design covariate files ----

NPCs <- 21
PCs <- paste0("PC", 1:NPCs)

# Get weight and BMI data to calculate baseline BMI from 
# weight or BMI measurement closest to baseline event date
QCd_for_baseline <- adiposity$weight
QCd_for_baseline$type <- "weight"
tmp <- adiposity$BMI 
tmp$type <- "BMI"
QCd_for_baseline <- bind_rows(QCd_for_baseline, tmp)
QCd_for_baseline <- QCd_for_baseline %>% arrange(eid, event_dt)
QCd_for_baseline <- split(QCd_for_baseline[, c("event_dt", "type", "value")], 
                          f = QCd_for_baseline$eid)

df <- adiposity$WHR
calc_covars <- df %>% group_by(eid) %>% 
  arrange(eid, event_dt) %>%
  summarise(baseline_date = first(event_dt),
            baseline_age = first(age_event),
            age_sq = baseline_age^2,
            FUyrs = interval(first(event_dt), last(event_dt)) / years(1),
            FU_n = n(),
            baseline_trait = first(value))

# Get baseline weight or BMI (nearest) from baseline file
calc_covars$baseline_adipo <- sapply(1:dim(calc_covars)[1], function (i) {
  dat_sub <- QCd_for_baseline[[as.character(calc_covars$eid[i])]]
  if (!is.null(dat_sub)) {
    closest_measure <- which.min(abs(calc_covars$baseline_date[i] - 
                                       dat_sub$event_dt))[1]
    # If the value is weight, keep as is, flip BMI to negative so we know
    # not to calculate BMI again
    res <- ifelse(dat_sub$type[closest_measure] == "weight", 
                  dat_sub$value[closest_measure],
                  -dat_sub$value[closest_measure])
  } else res <- NA
  return (res)
})

# Merge with previously calculated covariates
cleaned <- merge(calc_covars, general_covars, by = "eid")
# Set missing and inconsistent ancestry to "other"
cleaned$ancestry <- ifelse(cleaned$ancestry == "missing or inconsistent",
                           "other", cleaned$ancestry)
cleaned <- subset(cleaned, !is.na(cleaned$height))

# Convert baseline weight to baseline BMI
# Calculate BMI from weight (kg) and height (cm) or simply flip BMI sign
cleaned$baseline_BMI <- ifelse(cleaned$baseline_adipo < 0, 
                               -cleaned$baseline_adipo,
                               cleaned$baseline_adipo / (cleaned$height/100)^2) 

# Remove individuals missing any other covariate
QCd_covars <- cleaned[, c("eid", "sex", "ancestry",
                   "baseline_age", "age_sq", 
                   "height", "baseline_BMI", "baseline_trait",
                   "FUyrs", "FU_n", PCs)]
QCd_covars <- QCd_covars[complete.cases(QCd_covars), ]

# Calculate raw slopes ----

adiposity <- adiposity$WHR
covars <- QCd_covars

# One last check to remove any individuals without multiple measures
keep_eids <- covars$eid[covars$FU_n > 1]
dat <- subset(adiposity, adiposity$eid %in% keep_eids)

# mixed effects model for adiposity on age
mixed_model <- lmer(value ~ age_event + (age_event | eid), 
                    data = dat, REML = F)
# produce a new variable for the random slope
rs <- ranef(mixed_model)$eid$age_event
# create a new variable for the fixed effect of age on adiposity
bl <- fixef(mixed_model)[2]
# sum the random and fixed effect to get individual slope (BLUP)
slope <- rs + bl
# create slopes to add to covariate summary data frame
slopes <- data.frame(eid = rownames(ranef(mixed_model)$eid),
                     raw_slope = slope)

# Flag outlier slopes > 5 S.D. away from mean
# Remove values +/- 5 S.D. away from the strata mean slope
popn_mean <- mean(slopes$raw_slope, na.rm = T)
popn_sd <- sd(slopes$raw_slope, na.rm = T)
max_outlier <- popn_mean + 5*popn_sd
min_outlier <- popn_mean - 5*popn_sd

slopes$slope_outlier_flag <- slopes$raw_slope > max_outlier |
  slopes$raw_slope < min_outlier

# Add to individual-level covariates data
rs <- merge(covars, slopes, by = "eid")
# Exclude outlier-flagged slopes
rs <- subset(rs, !rs$slope_outlier_flag)

# Plot mean trajectory by slope quartile -----

rs$q <- cut(rs$raw_slope, quantile(rs$raw_slope), include.lowest = T,
            labels = paste0("q", 1:4))
df <- adiposity
df <- subset(df, df$eid %in% rs$eid)
df$q <- rs$q[match(df$eid, rs$eid)]

# Calculate mean and SE in each 5-year interval within each quartile
age_bin_cuts <- seq(20, 80, by = 5)
df$age_bin <- cut(df$age_event, age_bin_cuts, include.lowest = T)
plot_df <- df %>% group_by(q, age_bin) %>% 
  summarise(count = n(),
            mean_value = mean(value),
            se_value = sd(value)/sqrt(count))

# Plot 
p <- ggplot(plot_df, aes(x = age_bin, y = mean_value,
                         group = q, color = q, fill = q)) +
  geom_point() +
  geom_path() +
  geom_ribbon(aes(ymin = mean_value - se_value, 
                  ymax = mean_value + se_value),
              alpha = 0.2) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity")

pdf("tmp.pdf", onefile = T)
p
dev.off()

