# Author: Samvida S. Venkatesh
# Date: 11/10/21

library(tidyverse)

# Read data ----

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

resdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2211_models/traits_for_gwas/"
dir.create(resdir)

# BLUP files
blups <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SEX_STRATA, function (sx) {
    df <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2211_models/lmm_models/", 
                            p, "_", sx, "_blups_full_model.txt"),
                     sep = "\t", header = T, stringsAsFactors = F)
    colnames(df) <- c("eid", "b0", "b1")
    return (df)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(blups) <- PHENOTYPES

# IDs that passed sample QC
ids_passed_qc <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SEX_STRATA, function (sx) {
    read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/sample_qc/", 
                      p, "_", sx, "_ids_passed_qc.txt"),
               sep = "\t", header = T)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(ids_passed_qc) <- PHENOTYPES

# Covariates files (trait-specific) 
covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[PHENOTYPES]
general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

COVARS_LIST <- c("baseline_age", "age_sq", "year_of_birth", "FUyrs", "FU_n", 
                 "UKB_assmt_centre")

# Extract BLUP terms and relevant covariates ----

full_dat <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    # Filter to samples that passed QC (previously calculated)
    res <- blups[[p]][[sx]] %>% 
      filter(eid %in% ids_passed_qc[[p]][[sx]]$IID)
    # Merge in covariates
    res <- merge(res, covars[[p]][, c("eid", "baseline_age", "age_sq", 
                                      "FUyrs", "FU_n")],
                 by = "eid")
    # Merge in sex and UKB assessment centre
    res <- merge(res, 
                 general_covars[, c("eid", "sex", "year_of_birth", 
                                    "UKB_assmt_centre")], 
                 by = "eid")
    return (res)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(full_dat) <- PHENOTYPES

# Adjust coefficients for covariates and RINT ----

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    
    df <- full_dat[[p]][[sx]]
    # For sex-combined strata, adjust for sex
    adj_list <- COVARS_LIST
    if (sx == "sex_comb") adj_list <- c(adj_list, "sex")
    
    per_term <- lapply(c("b0", "b1"), function (blup_term) {
      
      # For slope, further adjust for intercept
      if (blup_term == "b1") adj_list <- c(adj_list, "b0")
      
      # Formula for adjustment
      mod_formula <- 
        formula(paste0(blup_term, " ~ ", paste(adj_list, collapse = " + ")))
      # Get residuals
      mod_resid <- lm(mod_formula, data = df)$residuals
      # RINT
      rinted_resid <- qnorm((rank(mod_resid) - 0.5) / sum(!is.na(mod_resid)))
      # Return formatted data
      res <- data.frame(FID = df$eid, 
                        IID = df$eid,
                        adj_trait = rinted_resid) 
      
      # Write results to table
      write.table(res,
                  paste0(resdir, p, "_", sx, "_", blup_term, ".txt"),
                  sep = "\t", row.names = F, quote = F)
    })
  })
})

