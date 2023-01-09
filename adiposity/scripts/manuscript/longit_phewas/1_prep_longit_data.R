# Author: Samvida S. Venkatesh
# Date: 31/08/21

library(lme4)
library(tidyverse)
theme_set(theme_bw())

# Read data ----

infile_path <- "" # REDACTED
gen_resources_path <- "" # REDACTED
outfile_path <- "" # REDACTED

PHENOS <- read.table(paste0(infile_path, "/code_lists/qcd_traits_available.txt"),
                     sep = "\t", header = F, stringsAsFactors = F)$V1
# Remvoe adiposity phenotypes
REMOVE_PHENOS <- c("BMI", "WC", "Weight", "WHR")

PHENOS <- PHENOS[!PHENOS %in% REMOVE_PHENOS]
SEX_STRATA <- c("F", "M", "sex_comb")
# Add sex as covariate for sex-combined analyses
MOD_COVARS <- c("baseline_age", "age_sq")
# Add data provider as covariate if there is more than one data provider
ADD_COVARS <- c("year_of_birth")

dat <- readRDS(paste0(infile_path, "/data/indiv_qcd_data.rds"))[PHENOS]
covars <- readRDS(paste0(infile_path, "/data/covariates.rds"))[PHENOS]

general_covars <- read.table(paste0(gen_resources_path, "/220504_QCd_demographic_covariates.txt"),
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

# Wrangle data ----

# Add in covariates
model_dat <- lapply(PHENOS, function (p) {
  
  res <- dat[[p]] %>% left_join(covars[[p]]) %>%
    mutate(t = age_event - baseline_age) %>%
    left_join(general_covars[, c("eid", "sex", "year_of_birth")],
                   by = "eid")
  
  return (res)
})
names(model_dat) <- PHENOS

# Create individual-level slopes ----

makeFormula <- function (dat_to_model) {
  # Covariates 
  include_covars <- c(MOD_COVARS, ADD_COVARS)
  # Include sex as covariate?
  if (length(unique(dat_to_model$sex)) > 1) 
    include_covars <- c(include_covars, "sex")
  # Include data provider as covariate?
  if (length(unique(dat_to_model$data_provider)) > 1)
    include_covars <- c(include_covars, "data_provider")
  
  # Write model formula 
  # With fixed effect of time and random effect of time
  mod_form <- paste0("value ~ t + ",
                     paste0(include_covars, collapse = " + "),
                     " + (t | eid)")
  return (mod_form)
}

full_models <- lapply(PHENOS, function (p) {
  cat(paste0("Running phenotype: ", p, "\n"))
  res <- lapply(SEX_STRATA, function (sx) {
    cat(paste0("\t", "in sex strata: ", sx, "\n"))
    sub_dat <- model_dat[[p]]
    if (sx != "sex_comb") 
      sub_dat <- model_dat[[p]] %>% filter(sex == sx)
    
    # Get formula
    mod_form <- makeFormula(sub_dat)
    
    # Run models
    lmod <- lmer(formula(mod_form), data = sub_dat, REML = F)
    return (lmod)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(full_models) <- PHENOS

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

lapply(PHENOS, function (p) {
  lapply(SEX_STRATA, function (sx) {
    full_blups <- getBLUPS(full_models[[p]][[sx]])
    write.table(full_blups, 
                paste0(outfile_path, "/longit_phewas/lmm_models/",
                       p, "_", sx, "_blups_full_model.txt"),
                sep = "\t", row.names = F, col.names = T, quote = F)
  })
})
