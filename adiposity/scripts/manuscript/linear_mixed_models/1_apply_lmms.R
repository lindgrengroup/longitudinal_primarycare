# Author: Samvida S. Venkatesh
# Date: 05/10/21

library(lme4)
library(tidyverse)
theme_set(theme_bw())

# Read files ----

mainpath <- "" # REDACTED
gen_resources_path <- "" # REDACTED
lmm_mods_path <- "" # REDACTED

# Parse in phenotype argument
args <- commandArgs(trailingOnly = T)
PHENO <- args[1]
SEX_STRATA <- c("F", "M", "sex_comb")

# Add sex as covariate for sex-combined analyses
MOD_COVARS <- c("baseline_age", "age_sq")
# Add data provider as covariate if there is more than one data provider
ADD_COVARS <- c("year_of_birth")

dat <- readRDS(paste0(mainpath, "/indiv_qcd_data.rds"))[[PHENO]]
dat$eid <- as.character(dat$eid)
covars <- readRDS(paste0(mainpath, "/covariates.rds"))[[PHENO]]
covars$eid <- as.character(covars$eid)

general_covars <- read.table(paste0(gen_resources_path, "/220504_QCd_demographic_covariates.txt"),
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

# Wrangle data ----

# Get time from baseline measurement and add in covariates
model_dat <- left_join(dat, 
                   covars[, c("eid", "baseline_age", "age_sq")], 
                   by = "eid")
model_dat <- left_join(model_dat, 
                       general_covars[, c("eid", 
                                          "sex", "year_of_birth")],
                       by = "eid")
model_dat <- model_dat %>% 
  mutate(t = age_event - baseline_age)

# Create individual-level slopes ----

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

minimal_models <- lapply(SEX_STRATA, function (sx) {
  sub_dat <- model_dat
  if (sx != "sex_comb") 
    sub_dat <- model_dat %>% filter(sex == sx)
  
  # Get formula
  mod_form <- makeFormula(sub_dat, add_adjust = F)
  
  # Run models
  lmod <- lmer(formula(mod_form), data = sub_dat, REML = F)
  return (lmod)
})
names(minimal_models) <- SEX_STRATA

# Save minimal_models
saveRDS(minimal_models, 
        paste0(lmm_mods_path, PHENO, "_minimal_model.rds"))

full_models <- lapply(SEX_STRATA, function (sx) {
  sub_dat <- model_dat
  if (sx != "sex_comb") 
    sub_dat <- model_dat %>% filter(sex == sx)
  
  # Get formula
  mod_form <- makeFormula(sub_dat, add_adjust = T)
  
  # Run models
  lmod <- lmer(formula(mod_form), data = sub_dat, REML = F)
  return (lmod)
})
names(full_models) <- SEX_STRATA

# Save minimal_models
saveRDS(full_models, 
        paste0(lmm_mods_path, PHENO, "_full_model.rds"))

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

model_blups <- lapply(SEX_STRATA, function (sx) {
  min_blups <- getBLUPS(minimal_models[[sx]])
  write.table(min_blups, 
              paste0(lmm_mods_path, PHENO, "_", sx, "_blups_minimal_model.txt"),
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  full_blups <- getBLUPS(full_models[[sx]])
  write.table(full_blups, 
              paste0(lmm_mods_path, PHENO, "_", sx, "_blups_full_model.txt"),
              sep = "\t", row.names = F, col.names = T, quote = F)
})
