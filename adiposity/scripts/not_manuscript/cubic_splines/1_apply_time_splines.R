# Author: Samvida S. Venkatesh
# Date: 03/11/21

library(lme4)
library(splines)
library(tidyverse)
theme_set(theme_bw())

# Parse in phenotype argument and analysis type - full or minimal
args <- commandArgs(trailingOnly = T)
PHENO <- args[1]
ATYPE <- args[2] 
SEX_STRATA <- c("F", "M", "sex_comb")

# Add sex as covariate for sex-combined analyses
MOD_COVARS <- c("baseline_age", "age_sq")
# Add data provider as covariate if there is more than one data provider
ADD_COVARS <- c("year_of_birth", "smoking_status")

dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/indiv_qcd_data.rds")[[PHENO]]
dat$eid <- as.character(dat$eid)
covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[[PHENO]]
covars$eid <- as.character(covars$eid)

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220131_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

# Wrangle data ----

# Get time from baseline measurement and add in covariates
model_dat <- left_join(dat, 
                       covars[, c("eid", "baseline_age", "age_sq")], 
                       by = "eid")
model_dat <- left_join(model_dat, 
                       general_covars[, c("eid", 
                                          "sex", "year_of_birth", "smoking_status")],
                       by = "eid")
model_dat <- model_dat %>% 
  mutate(t = age_event - baseline_age)

# Functions to run natural cubic spline models ---- 

# Run models with 3:10 df in fixed effect and 3 df in random effect 
# to choose best models
# If all the models are NA that's because we have too many random effects
# but not enough observations, so keep reducing the df of random effects

# Function to create model formula
getModFormula <- function (dat_to_model, ndf_fe, ndf_re, add_adjust = F) {
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
  
  mod_formula <- formula(paste0("value ~ ", 
                                paste0(include_covars, collapse = " + "), 
                                "+ ns(t, df = ", ndf_fe, 
                                ") + (ns(t, df = ", ndf_re, 
                                ") | eid)"))
  return (mod_formula)
}

# Function to run through all possible degrees of freedom for 
# given set of data
runThroughMods <- function (dat_to_model, add_adjust = F, 
                            min_fe = 3, max_fe = 10, 
                            max_re = 3, 
                            log_file) {
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
      res <- tryCatch(lmer(getModFormula(dat_to_model, 
                                         ndf_fe, ndf_re, 
                                         add_adjust = add_adjust),
                           data = dat_to_model, REML = F),
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

models <- lapply(SEX_STRATA, function (sx) {
  sub_dat <- model_dat
  if (sx != "sex_comb") 
    sub_dat <- model_dat %>% filter(sex == sx)
  
  # Log file
  log_file <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/spline_models/logs/",
                     PHENO, "_", sx, "_", ATYPE, "_models.txt")
  # This will error out if ATYPE is wrong
  add_adjustments <- ifelse(ATYPE == "full", T, 
                            ifelse(ATYPE == "minimal", F,
                                                       NA))
  mod_list <- runThroughMods(sub_dat, 
                             add_adjust = add_adjustments,
                             min_fe = 3, max_fe = 10, 
                             max_re = 3, 
                             log_file)
  bics <- sapply(mod_list, function (x) BIC(x))
  best_mod <- mod_list[[which.min(bics)]]
  
  # Log results
  sink(log_file, append = T) 
  cat(paste0("Best ", ATYPE, " model in strata ", sx, ", phenotype: ", PHENO, "\n"))
  print(summary(best_mod))
  cat(paste0("########################################################", "\n"))
  sink()
  
  # Save model
  return (best_mod)
})
names(models) <- SEX_STRATA

# Save models
saveRDS(models, 
        paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/spline_models/",
               PHENO, "_", ATYPE, "_model.rds"))

# Save coefficient values (fixef + ranef) for individuals ----

getBLUPS <- function (mod) {
  res <- coef(mod)$eid
  # We only want the columns with both fixed and random effects 
  # (intercept, "t", and any cubic spline terms of t)
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
  blups <- getBLUPS(models[[sx]])
  write.table(blups, 
              paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/spline_models/",
                     PHENO, "_", sx, "_blups_", ATYPE, "_model.txt"),
              sep = "\t", row.names = F, col.names = T, quote = F)
})
