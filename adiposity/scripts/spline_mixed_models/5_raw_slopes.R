# Author: Samvida S. Venkatesh
# Date: 13/05/2021

library(lme4)
library(gamm4)
library(tidyverse)

# Read files ----

adiposity <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/visually_QCd_adiposity.rds")
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/model_covariates.rds")

# Read arguments and prepare subsets for analysis ----

args <- commandArgs(trailingOnly = T)

PHENO <- args[1]
STRATUM <- args[2]
STRATUM_ANCESTRY <- sub("_.*", "", STRATUM)
STRATUM_SEX <- sub(".*_", "", STRATUM)

sub_covars <- covars[[PHENO]]
adipo <- adiposity[[PHENO]]

# Subset adiposity data relevant to stratum
# One last check to remove any individuals without multiple measures
si_eids <- sub_covars$eid[sub_covars$sex == STRATUM_SEX &
                            sub_covars$ancestry == STRATUM_ANCESTRY &
                            sub_covars$FU_n > 1]
dat <- subset(adipo, adipo$eid %in% si_eids)

# Create individual-level raw slopes ----

## linear mixed effects model ----

linear_model <- lmer(value ~ age_event + (1 + age_event | eid),
                     data = dat)

# Save linear model
saveRDS(linear_model, 
        paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/linear/", PHENO,
               "_linear_model_", STRATUM, ".rds"))


## mixed effects model with splines ----

NKNOTS <- ifelse(PHENO == "BMI", 8, 3)
spline_model <- gamm4(value ~ s(age_event, k = NKNOTS, bs = "cr"), 
                      random = ~(1 + age_event | eid), 
                      data = dat)

# Save spline model
saveRDS(spline_model, 
        paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/splines/", PHENO,
               "_spline_model_", STRATUM, ".rds"))

## create covariates and slopes data-frame ----

# produce a new variable for the random slope
rs <- ranef(spline_model$mer)$eid$age_event
# create a new variable for the fixed effect of age on adiposity
bl <- fixef(spline_model$mer)[2]
# sum the random and fixed effect to get individual slope (BLUP)
spline_slope <- rs + bl
# create slopes to add to covariate summary data frame
spline_slopes <- data.frame(eid = rownames(ranef(spline_model$mer)$eid),
                     spline_slope = spline_slope)

rs <- ranef(linear_model)$eid$age_event
bl <- fixef(linear_model)[2]
linear_slope <- rs + bl
# create slopes to add to covariate summary data frame
linear_slopes <- data.frame(eid = rownames(ranef(linear_model)$eid),
                     linear_slope = linear_slope)

# Merge slopes dataframe
slopes <- merge(linear_slopes, spline_slopes, by = "eid")

# Flag outlier slopes > 5 S.D. away from mean
popn_mean <- mean(slopes$linear_slope, na.rm = T)
popn_sd <- sd(slopes$linear_slope, na.rm = T)
max_outlier <- popn_mean + 5*popn_sd
min_outlier <- popn_mean - 5*popn_sd
slopes$linear_slope_outlier_flag <- slopes$linear_slope > max_outlier |
  slopes$linear_slope < min_outlier

popn_mean <- mean(slopes$spline_slope, na.rm = T)
popn_sd <- sd(slopes$spline_slope, na.rm = T)
max_outlier <- popn_mean + 5*popn_sd
min_outlier <- popn_mean - 5*popn_sd
slopes$spline_slope_outlier_flag <- slopes$spline_slope > max_outlier |
  slopes$spline_slope < min_outlier

# Report QC metrics
sink(paste0("log_files/splines/spline_slopes_QC_", PHENO, ".txt"), append = T)
cat(paste0("Strata: ", 
           STRATUM_ANCESTRY, " ", STRATUM_SEX, "\n",
           "**FILTER** EXCLUDED, Raw slope > 5 S.D. away from strata mean: ", 
           sum(slopes$spline_slope_outlier_flag), "\n"))
sink()

# Add to individual-level covariates data
slopes$eid <- as.character(slopes$eid)
res <- merge(sub_covars, slopes, by = "eid")

# Save covariates with slopes
saveRDS(res, 
        paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/", PHENO,
               "_raw_slopes_and_covars_", STRATUM, ".rds"))