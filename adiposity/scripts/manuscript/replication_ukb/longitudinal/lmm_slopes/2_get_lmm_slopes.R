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

dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/data/main_data_adipo_change.rds")[PHENO]

discovery_indivs <- lapply(PHENO, function (p) {
  res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/traits_for_GWAS/lmm_traits_", 
                           p, "_sex_comb.txt"), 
                    sep = "\t", header = T, stringsAsFactors = F)$IID
  return (res)
})
names(discovery_indivs) <- PHENO

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

indivs_dementia <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eids_with_dementia.txt",
                              sep = "\t", header = F, stringsAsFactors = F)$V1
indivs_dementia <- as.character(indivs_dementia)

# Wrangle data ----

# Exclude individuals in discovery strata
# Get time from baseline measurement and add in covariates
model_dat <- lapply(PHENO, function (p) {
  res <- dat[[p]] %>% filter(!eid %in% discovery_indivs[[p]])
  
  res <- res %>%
    group_by(eid) %>% arrange(age_event, .by_group = T) %>%
    mutate(baseline_age = first(age_event),
           age_sq = baseline_age^2,
           t = age_event - baseline_age)
  res <- left_join(res, general_covars[, c("eid", 
                                           "sex", "year_of_birth")],
                   by = "eid")
  res$dementia_status <- res$eid %in% indivs_dementia
  return (res)
})
names(model_dat) <- PHENO

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

full_models <- lapply(PHENO, function (p) {
  res <- lapply(SEX_STRATA, function (sx) {
    sub_dat <- model_dat[[p]]
    if (sx != "sex_comb") 
      sub_dat <- model_dat[[p]] %>% filter(sex == sx)
    # Get subset of data without dementia
    sub_dat_no_dementia <- sub_dat %>% filter(!dementia_status)
    
    # Get formula
    mod_form <- makeFormula(sub_dat, add_adjust = T)
    
    # Run models
    lmod_all <- lmer(formula(mod_form), data = sub_dat, REML = F)
    lmod_no_dementia <- lmer(formula(mod_form), data = sub_dat_no_dementia, 
                             REML = F)
    
    return (list(lmod_all = lmod_all,
                 lmod_no_dementia = lmod_no_dementia))
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(full_models) <- PHENO
# Save models
saveRDS(full_models, 
        "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/lmm_models/full_models.rds")

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
  lapply(SEX_STRATA, function (sx) {
    all_blups <- getBLUPS(full_models[[p]][[sx]]$lmod_all)
    write.table(all_blups, 
                paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/lmm_models/",
                       p, "_", sx, "_all_blups.txt"),
                sep = "\t", row.names = F, col.names = T, quote = F)
    
    blups_no_dementia <- getBLUPS(full_models[[p]][[sx]]$lmod_no_dementia)
    write.table(blups_no_dementia, 
                paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/lmm_models/",
                       p, "_", sx, "_blups_no_dementia.txt"),
                sep = "\t", row.names = F, col.names = T, quote = F)
  })
})
