# Author: Samvida S. Venkatesh
# Date: 05/10/21

library(lme4)
library(tidyverse)
theme_set(theme_bw())

# Read files ----

mainpath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data"
gpdat_path <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity"
gen_resources_path <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources"
outpath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sensitivity_rs429358"

# Parse in phenotype argument
args = commandArgs(trailingOnly=TRUE)
PHENO <- args[1]
SEX_STRATA <- c("F", "M", "sex_comb")

# Add sex as covariate for sex-combined analyses
MOD_COVARS <- c("baseline_age", "age_sq", "year_of_birth", "data_provider",
                paste0("PC", 1:21))
NO_AGE_COVARS <- c("year_of_birth", "data_provider",
                   paste0("PC", 1:21))

dat <- readRDS(paste0(mainpath, "/indiv_qcd_data.rds"))[[PHENO]]
dat$eid <- as.character(dat$eid)
covars <- readRDS(paste0(mainpath, "/covariates.rds"))[[PHENO]]
covars$eid <- as.character(covars$eid)

general_covars <- read.table(paste0(gen_resources_path, "/220504_QCd_demographic_covariates.txt"),
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

# Genotypes / dosages at rs429358
apoe_dosages <- read.table(paste0(gpdat_path, "/sample_variant_counts/rs429358_dosages.txt"),
                           sep = " ", header = T, stringsAsFactors = F)
# Remove first row, which contains info on type of column and columns 
# 2, 3, 4 (ID repeat, missingness, sex)
apoe_dosages <- apoe_dosages[-1, c(1, 5)]
colnames(apoe_dosages) <- c("eid", "dosage")
apoe_dosages$eid <- as.character(apoe_dosages$eid)

# Wrangle data ----

# Get time from baseline measurement and add in covariates
model_dat <- left_join(dat, 
                       covars[, c("eid", "baseline_age", "age_sq")], 
                       by = "eid")
model_dat <- left_join(model_dat, 
                       general_covars[, c("eid", 
                                          "sex", "year_of_birth", paste0("PC", 1:21))],
                       by = "eid")
model_dat <- model_dat %>% 
  mutate(t = age_event - baseline_age)

# Add in APOE genotype dosage
model_dat <- model_dat %>%
  left_join(apoe_dosages, by = "eid") %>%
  mutate(dosage = as.numeric(dosage))

model_dat <- model_dat[complete.cases(model_dat), ]

# Run models ----

full_models <- lapply(SEX_STRATA, function (sx) {
  sub_dat <- model_dat
  include_covs <- MOD_COVARS
  include_covs_no_age <- NO_AGE_COVARS
  if (sx != "sex_comb") {
    sub_dat <- model_dat %>% filter(sex == sx)
  } else {
    include_covs <- c(MOD_COVARS, "sex")
    include_covs_no_age <- c(NO_AGE_COVARS, "sex")
  }
  
  # Adjusted for age
  age_mod <- lmer(formula(paste0("value ~ t + dosage + dosage:t + ",
                                 paste0(include_covs, collapse = " + "),
                                 " + (t | eid)")), 
                  data = sub_dat)
  
  # Unadjusted for age
  no_age_mod <- lmer(formula(paste0("value ~ t + dosage + dosage:t + ",
                                    paste0(include_covs_no_age, collapse = " + "),
                                    " + (t | eid)")), 
                     data = sub_dat)
  
  # Get results to write in table
  age_mod_summ <- as.data.frame(summary(age_mod)$coefficients[c("dosage", "t:dosage"), 
                                                              c(1:3)])
  age_mod_summ$term <- c("SNP", "SNP:time")
  age_mod_summ$adjustment <- "age_adj"
  
  no_age_mod_summ <- as.data.frame(summary(no_age_mod)$coefficients[c("dosage", "t:dosage"), 
                                                              c(1:3)])
  no_age_mod_summ$term <- c("SNP", "SNP:time")
  no_age_mod_summ$adjustment <- "no_age_adj"
  
  res <- bind_rows(age_mod_summ, no_age_mod_summ)
  colnames(res) <- c("beta", "stderr", "tstat", "adjustment")
  res$sex_strata <- sx
  res$phenotype <- PHENO
  return (list(age_adj_mod = age_mod,
               no_age_adj_mod = no_age_mod,
               res_table = res))
})
names(full_models) <- SEX_STRATA

saveRDS(full_models, 
        paste0(outpath, "/full_models_for_", PHENO, ".rds"))

restable <- lapply(full_models, function (mlist) {
  return (mlist$res_table)
})
restable <- bind_rows(restable)

write.table(restable,
            paste0(outpath, "/lme_results_", PHENO, ".txt"),
            sep = "\t", quote = F, row.names = F)
