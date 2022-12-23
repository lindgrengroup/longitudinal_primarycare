# Author: Samvida S. Venkatesh
# Date: 22/12/2022

library(tidyverse)
library(broom)
theme_set(theme_bw())

RANDOM_SEED <- 221222
set.seed(RANDOM_SEED)

# Read files ----

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

# Adiposity data
dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/indiv_qcd_data.rds")[PHENOTYPES]
dat$eid <- as.character(dat$eid)

# Covariates
general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

# Lead SNPs from GWAS
vars_to_replicate <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/data/b0_snps_to_replicate.txt",
                                sep = "\t", header = T, stringsAsFactors = F)
VARIDS <- vars_to_replicate$SNP

# Genotypes / dosages at variants of interest
var_dosages <- lapply(VARIDS, function (varid) {
  res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sample_variant_counts/",
                           varid, "_dosages.txt"),
                    sep = " ", header = T, stringsAsFactors = F)
  # Remove first row, which contains info on type of column and columns 
  # 2, 3, 4 (ID repeat, missingness, sex)
  res <- res[-1, c(1, 5)]
  colnames(res) <- c("eid", varid)
  return (res)
})
names(var_dosages) <- VARIDS

# Samples that passed genotyping QC
# IDs that passed sample QC
ids_passed_qc <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SEX_STRATA, function (sx) {
    df <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/sample_qc/", 
                            p, "_", sx, "_ids_passed_qc.txt"),
                     sep = "\t", header = T)
    df$eid <- as.character(df$IID)
    return (df)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(ids_passed_qc) <- PHENOTYPES

MOD_COVARS <- c("sex", "age_event", "age_sq", 
                "year_of_birth", "data_provider", "UKB_assmt_centre",
                paste0("PC", 1:21), "genotyping.array") 
general_covars <- general_covars %>%
  dplyr::select(any_of(c("eid", MOD_COVARS)))

# Wrangle data ----

# Convert dosages to 0/1/2 genotypes based on threshold
DOSAGE_THRESHOLD_0 <- 0.5
DOSAGE_THRESHOLD_2 <- 1.5
var_dosages_hardcall <- lapply(VARIDS, function (varid) {
  res <- var_dosages[[varid]] %>% 
    mutate(genotype = ifelse(!!as.symbol(varid) < DOSAGE_THRESHOLD_0, "0",
                             ifelse(!!as.symbol(varid) > DOSAGE_THRESHOLD_2, 
                                    "2", "1")),
           eid = as.character(eid)) %>%
    dplyr::select(all_of(c("eid", "genotype")))
  return (res)
})
names(var_dosages_hardcall) <- VARIDS

# Function to add hard-called genotype
addGenoGroup <- function (df, varid) {
  res <- df
  var_dat <- var_dosages_hardcall[[varid]]
  res$genotype <- var_dat$genotype[match(res$eid, var_dat$eid)]
  res <- res[!is.na(res$genotype), ]
  return (res)
}

# Testing functions ----

TRAIT_COVARS <- c("sex", "age_event", "age_sq", 
                  "year_of_birth", "data_provider", "UKB_assmt_centre") 
GENETIC_COVARS <- c(paste0("PC", 1:21), "genotyping.array")

quantTest <- function (df, ss) {
  covars_include <- TRAIT_COVARS
  if (ss != "sex_comb") {
    covars_include <- covars_include[-which(covars_include == "sex")]
  }
  if (length(unique(df$data_provider)) < 2)
    covars_include <- covars_include[-which(covars_include == "data_provider")]
  
  # Rank-based inverse normally transform the residuals of the trait
  # Formula for adjustment
  mod_formula <- 
    formula(paste0("value ~ ", paste0(covars_include, collapse = " + ")))
  # Get residuals
  mod_resid <- lm(mod_formula, data = df)$residuals
  # RINT
  rinted_resid <- qnorm((rank(mod_resid) - 0.5) / sum(!is.na(mod_resid)))
  # Return formatted data
  resid_dat <- data.frame(eid = df$eid, 
                    adj_trait = rinted_resid) 
  stage2_df <- df %>% inner_join(resid_dat, by = "eid")
  
  if (ss != "sex_comb") {
    mod_formula <- 
      formula(paste0("adj_trait ~ genotype + ", paste0(GENETIC_COVARS, collapse = " + ")))
  } else {
    mod_formula <- 
      formula(paste0("adj_trait ~ genotype + sex + ", paste0(GENETIC_COVARS, collapse = " + ")))
  }
  
  modeled_dat <- lm(mod_formula, data = stage2_df)
  
  print_res <- tidy(modeled_dat) %>%
    filter(term == "genotype") %>%
    rename(beta = estimate, se = std.error, tstat = statistic, pval = p.value) %>%
    dplyr::select(all_of(c("beta", "se", "tstat", "pval")))
  
  return (print_res)
}

# Apply to each strata and variant ----

FAC_COVARS <- c("sex", "data_provider", "UKB_assmt_centre", "genotyping.array")
NUM_COVARS <- c("age_event", "age_sq", "year_of_birth", 
                paste0("PC", 1:21))

full_res <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    # Join data with covariates and calculate additional covariates
    df <- dat[[p]] %>% left_join(general_covars, by = "eid") 
    df <- df %>% inner_join(ids_passed_qc[[p]][[sx]])
    
    # Get only baseline value for each individual
    df <- df %>% group_by(eid) %>%
      arrange(age_event) %>% slice(1)
    
    df <- df %>%
      mutate(age_sq = age_event^2) %>%
      mutate(across(all_of(FAC_COVARS), factor)) %>%
      mutate(across(all_of(NUM_COVARS), as.numeric)) %>%
      dplyr::select(any_of(c("eid", "value", MOD_COVARS)))
    
    # Go through each variant to replicate
    per_varid <- lapply(VARIDS, function (v) {
      
      full_df <- addGenoGroup(df, v) %>%
        mutate(genotype = as.numeric(genotype))
      test_res <- quantTest(full_df, sx)
      test_res$varid <- v
      
      return (test_res)
    })
    per_varid <- bind_rows(per_varid)
    per_varid$sex_strata <- sx
    return (per_varid)
  })
  per_sex <- bind_rows(per_sex)
  per_sex$phenotype <- p
  return (per_sex)
})
full_res <- bind_rows(full_res)

write.table(full_res,
            "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2211_models/GWAS/post_GWAS/baseline_in_same_data_replication.txt",
            sep = "\t", row.names = F, quote = F)
