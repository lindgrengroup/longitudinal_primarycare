# Author: Samvida S. Venkatesh
# Date: 05/10/21

library(lme4)
library(tidyverse)
theme_set(theme_bw())

# Read files ----

# Parse in phenotype arguments
PHENO <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

# Add sex as covariate for sex-combined analyses
MOD_COVARS <- c("baseline_age", "age_sq")
# Add data provider as covariate if there is more than one data provider
ADD_COVARS <- c("year_of_birth")

dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/non_wb_gp_main_data_passed_longit_filter.rds")[PHENO]

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

ANCESTRIES <- unique(general_covars$ancestry)

# Wrangle data ----

# Get time from baseline measurement and add in covariates
model_dat <- lapply(PHENO, function (p) {
  
  res <- dat[[p]] %>%
    group_by(eid) %>% arrange(age_event, .by_group = T) %>%
    mutate(baseline_age = first(age_event),
           age_sq = baseline_age^2,
           t = age_event - baseline_age)
  res <- left_join(res, general_covars[, c("eid", 
                                           "sex", "year_of_birth", "ancestry")])
  return (res)
})
names(model_dat) <- PHENO

# Create individual-level slopes within each ancestry group and sex ----

makeFormula <- function (dat_to_model, add_adjust = F) {
  # Covariates 
  include_covars <- MOD_COVARS
  # Additional covariates?
  if (add_adjust) include_covars <- c(MOD_COVARS, ADD_COVARS)
  # Include sex as covariate?
  if (length(unique(dat_to_model$sex)) > 1) 
    include_covars <- c(include_covars, "sex")
  # Include data provider as covariate?
  # Only if add adjustments is also true
  if (length(unique(dat_to_model$data_provider)) > 1 & add_adjust)
    include_covars <- c(include_covars, "data_provider")
  
  # Write model formula 
  # With fixed effect of time and random effect of time
  mod_form <- paste0("value ~ t + ",
                     paste0(include_covars, collapse = " + "),
                     " + (t | eid)")
  return (mod_form)
}

full_models <- lapply(PHENO, function (p) {
  res_list <- lapply(ANCESTRIES, function (anc) {
    # Subset by ancestry
    sub_dat <- model_dat[[p]] %>% filter(ancestry == anc)
    res <- lapply(SEX_STRATA, function (sx) {
      # Subset by sex
      to_model <- sub_dat
      if (sx != "sex_comb") 
        to_model <- to_model %>% filter(sex == sx)

      # Get formula
      mod_form <- makeFormula(to_model, add_adjust = T)
      
      # Run model
      lmod <- lmer(formula(mod_form), data = to_model, REML = F)

      return (lmod)
    })
    names(res) <- SEX_STRATA
    return (res)
  })
  names(res_list) <- ANCESTRIES
  return (res_list)
})
names(full_models) <- PHENO
# Save models
saveRDS(full_models, 
        "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/non_wb_ancestry/lmm_models/full_models.rds")

# Save coefficient values (fixef + ranef) for individuals ----

getBLUPS <- function (mod) {
  res <- coef(mod)$eid
  # We only want the columns with both fixed and random effects 
  # (intercept and "t")
  # But we want to get the full effect (fixed + random) for it 
  
  # Remove all possible covariate columns
  to_remove <- c(MOD_COVARS, ADD_COVARS, "sex", "data_provider")
  remove_cols <- grep(paste0(to_remove, collapse = "|"),
                      colnames(res))
  res <- res[, -remove_cols]
  
  # Remove any remaining columns that only have fixed effects
  # because their variance will be 0
  keep_cols <- apply(res, 2, function (x) var(x) != 0)
  res <- res[, keep_cols]
  vars_kept <- colnames(res)
  
  # Return dataframe to write
  res$eid <- rownames(res)
  res <- res[, c("eid", vars_kept)]
  
  return (res)
}

lapply(PHENO, function (p) {
  lapply(ANCESTRIES, function (anc) {
    lapply(SEX_STRATA, function (sx) {
      all_blups <- getBLUPS(full_models[[p]][[anc]][[sx]])
      write.table(all_blups, 
                  paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/non_wb_ancestry/lmm_models/",
                         p, "_", anc, "_", sx, "_all_blups.txt"),
                  sep = "\t", row.names = F, col.names = T, quote = F)
    })
  })
})
