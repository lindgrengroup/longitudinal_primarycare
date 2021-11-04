# Author: Samvida S. Venkatesh
# Date: 05/10/21

library(lme4)
library(tidyverse)
theme_set(theme_bw())

# Read files ----

# Parse in phenotype argument
args <- commandArgs(trailingOnly = T)
PHENO <- args[1]
SEX_STRATA <- c("F", "M", "sex_comb")
PCs <- paste0("PC", 1:21)
COVARS <- c("baseline_age", "age_sq", "data_provider")

dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/indiv_qcd_data.rds")[[PHENO]]
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/covariates.rds")[[PHENO]]

# Wrangle data ----

# Get time from baseline measurement and add in covariates
model_dat <- merge(dat, covars, by = "eid")
model_dat <- model_dat %>% mutate(t = age_event - baseline_age)

# Create individual-level slopes ----

slope_models <- lapply(SEX_STRATA, function (sx) {
  if (sx != "sex_comb") {
    sub_dat <- model_dat %>% filter(sex == sx)
    mod_covars <- COVARS
  } else {
    sub_dat <- model_dat
    mod_covars <- c(COVARS, "sex")
  }
  # One last check to remove any individuals without multiple measures
  sub_dat <- sub_dat %>% filter(FU_n > 1)
  
  # mixed effects model for adiposity on age
  covars_form <- paste0("value ~ ", 
                        paste0(mod_covars, collapse = " + "), " + ",
                        paste0(PCs, collapse = " + "))
  full_form <- paste0(covars_form, " + t + (t | eid)")
  lmod <- lmer(formula(full_form), data = sub_dat, REML = F)
  return (lmod)
})
names(slope_models) <- SEX_STRATA

# Save slopes models
saveRDS(slope_models, 
        paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/lmm_",
               PHENO, ".rds"))

# Save coefficient values (fixef + ranef) for individuals ----

blups <- lapply(SEX_STRATA, function (sx) {
  res <- coef(slope_models[[sx]])$eid
  # Only get the columns with random effect terms as well as f.e.
  all_covs <- paste0(c(COVARS, PCs, "sex"), collapse = "|")
  remove_cols <- grep(all_covs, colnames(res))
  res <- res[, -remove_cols]
  # Remove any columns with only fixed effects
  keep_cols <- apply(res, 2, function (x) var(x) != 0)
  res <- res[, keep_cols]
  return (res)
})
names(blups) <- SEX_STRATA
saveRDS(blups, 
        paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/lmm_blups_",
               PHENO, ".rds"))
