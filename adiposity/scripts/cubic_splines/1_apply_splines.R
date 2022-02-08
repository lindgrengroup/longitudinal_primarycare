# Author: Samvida S. Venkatesh
# Date: 03/11/21

library(lme4)
library(splines)
library(tidyverse)
theme_set(theme_bw())

# Parse in phenotype argument
args <- commandArgs(trailingOnly = T)
PHENO <- args[1]
SEX_STRATA <- c("F", "M", "sex_comb")
COVARS <- c("baseline_age", "age_sq")

dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/indiv_qcd_data.rds")[[PHENO]]
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/covariates.rds")[[PHENO]]

# Wrangle data ----

# Get time from baseline measurement and add in covariates
model_dat <- merge(dat, covars, by = "eid")
model_dat <- model_dat %>% mutate(t = age_event - baseline_age,
                                  age_event_sq = age_event^2)

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
                                paste0(mod_covars, collapse = " + "), 
                                "+ ns(t, df = ", ndf_fe, 
                                ") + (ns(t, df = ", ndf_re, 
                                ") | eid)"))
  return (mod_formula)
}

# Function to run through all possible dfs
runThroughMods <- function (sx, min_fe = 3, max_fe = 10, 
                            max_re = 3, log_file) {
  # Get data to model
  sub_dat <- model_dat
  if (sx != "sex_comb") sub_dat <- sub_dat %>% filter(sex == sx)
  
  ndf_re <- max_re
  all_mods <- list()
  # Step through models until at least one produces a non-NA result
  while (all(is.na(all_mods)) & ndf_re > 0) {
    
    sink(log_file, append = T) 
    print(paste0("Testing ndf_re: ", ndf_re))
    sink()
    
    # Step through all df for fixed effects
    all_mods <- lapply(min_fe:max_fe, function (ndf_fe) {
      
      sink(log_file, append = T)
      print(paste0("        Running ndf_fe: ", ndf_fe))
      sink()
      
      # Get formula for given dfs and apply
      res <- tryCatch(lmer(getModFormula(sx, ndf_fe, ndf_re),
                           data = sub_dat, REML = F),
                      error = function (err) NA)
      
      # Print BIC of model
      sink(log_file, append = T)
      if (!is.na(res)) { print(paste0("        BIC: ", BIC(res))) }
      else { print("        Model failed")}
      sink()
      
      return (res)
    })
    ndf_re <- ndf_re - 1
  }
  # Return models that are not NA
  all_mods <- all_mods[!is.na(all_mods)]
  return (all_mods)
}

# Apply to all sex strata and log results ----

spline_mods <- lapply(SEX_STRATA, function (sx) {
  # Log file
  log_file <- paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/log_files/age_adj_time_models_",
                     PHENO, "_", sx, ".txt")
  mod_list <- runThroughMods(sx, min_fe = 3, max_fe = 10, 
                             max_re = 3, log_file)
  bics <- sapply(mod_list, function (x) BIC(x))
  best_mod <- mod_list[[which.min(bics)]]
  
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
        paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/age_adj_",
               PHENO, ".rds"))

# Save coefficient values (fixef + ranef) for individuals ----

blups <- lapply(SEX_STRATA, function (sx) {
  res <- coef(spline_mods[[sx]])$eid
  # Only get the columns with random effect terms as well as f.e.
  all_covs <- paste0(c(COVARS, "sex"), collapse = "|")
  remove_cols <- grep(all_covs, colnames(res))
  res <- res[, -remove_cols]
  # Remove any columns with only fixed effects
  keep_cols <- apply(res, 2, function (x) var(x) != 0)
  res <- res[, keep_cols]
  return (res)
})
names(blups) <- SEX_STRATA
saveRDS(blups, 
        paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/age_adj_blups_",
               PHENO, ".rds"))
