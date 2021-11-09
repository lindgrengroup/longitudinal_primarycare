# Author: Samvida S. Venkatesh
# Date: 03/11/21

library(lme4)
library(splines)
library(tidyverse)

# Parse in phenotype argument
args <- commandArgs(trailingOnly = T)
PHENO <- args[1]
SEX_STRATA <- c("F", "M", "sex_comb")
PCs <- paste0("PC", 1:21)
COVARS <- c("data_provider")

dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/indiv_qcd_data.rds")[[PHENO]]
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/covariates.rds")[[PHENO]]

log_file <- paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/log_files/cubic_splines_",
                   PHENO, ".txt")

# Wrangle data ----

# Add in covariates
model_dat <- merge(dat, covars, by = "eid")

# Functions to run natural cubic spline models ---- 

# Run models with 3:10 df in fixed effect and 3 df in random effect 
# to choose best models
# If all the models are NA that's because we have too many random effects
# but not enough observations, so keep reducing the df of random effects

# Function to create model formula
getModFormula <- function (sx, ndf_fe, ndf_re) {
  # Add sex as covariate to sex-combined models
  mod_covars <- COVARS
  if (sx == "sex_comb") mod_covars <- c(mod_covars, "sex")
  
  mod_formula <- formula(paste0("value ~ ", 
                                paste0(mod_covars, collapse = " + "), " + ", 
                                paste0(PCs, collapse = " + "), 
                                "+ ns(age_event, df = ", ndf_fe, 
                                ") + (ns(age_event, df = ", ndf_re, 
                                ") | eid)"))
  return (mod_formula)
}

# Function to run through all possible dfs
runThroughMods <- function (sx, min_fe = 3, max_fe = 10, 
                            max_re = 3) {
  # Get data to model
  sub_dat <- model_dat
  if (sx != "sex_comb") sub_dat <- sub_dat %>% filter(sex == sx)
  
  ndf_re <- max_re
  all_mods <- list()
  # Step through models until at least one produces a non-NA result
  while (all(is.na(all_mods)) & ndf_re > 0) {
    print(paste0("Testing ndf_re: ", ndf_re))
    # Step through all df for fixed effects
    all_mods <- lapply(min_fe:max_fe, function (ndf_fe) {
      print(paste0("        Running ndf_fe: ", ndf_fe))
      # Get formula for given dfs and apply
      res <- tryCatch(lmer(getModFormula(sx, ndf_fe, ndf_re),
                           data = sub_dat, REML = F),
                      error = function (err) NA)
      return (res)
    })
    ndf_re <- ndf_re - 1
  }
  # Return models that are not NA
  all_mods <- all_mods[!is.na(all_mods)]
  return (all_mods)
}

# Function to pick model with lowest BIC (parsimonious)
pickBestMod <- function (mod_list) {
  bics <- sapply(mod_list, function (x) BIC(x))
  best_mod <- mod_list[[which.min(bics)]]
  return (best_mod)
}

# Apply to all sex strata and log results ----

spline_mods <- lapply(SEX_STRATA, function (sx) {
  mod_list <- runThroughMods(sx, min_fe = 3, max_fe = 10, 
                             max_re = 3)
  best_mod <- pickBestMod(mod_list)
  # Log results
  sink(log_file, append = T) 
  cat(paste0("Best model in strata ", sx, ", phenotype: ", PHENO, "\n"))
  print(summary(best_mod))
  cat(paste0("########################################################", "\n"))
  sink()
  # Save model
  return (best_mod)
})
names(spline_mods) <- SEX_STRATA

saveRDS(spline_mods, 
        paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines_",
               PHENO, ".rds"))

# Save coefficient values (fixef + ranef) for individuals ----

blups <- lapply(SEX_STRATA, function (sx) {
  res <- coef(spline_mods[[sx]])$eid
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
        paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_spline_blups_",
               PHENO, ".rds"))
