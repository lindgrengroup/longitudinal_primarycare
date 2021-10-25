# Author: Samvida S. Venkatesh
# Date: 11/10/21

library(tidyverse)

# Read data ----

PHENOTYPES <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/pheno_names.txt",
                         sep = "\t", header = F, stringsAsFactors = F)$V1

# BLUPs
blups <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/lmm_blups_all_strata.rds")
PHENOTYPES <- names(blups)
SEX_STRATA <- names(blups[[1]])

# IDs that passed sample QC for GWAS 
qcd_ids <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/ids_passed_qc_211015.txt", 
                      sep = "\t", header = T)

# Covariates file (general and trait-specific) 
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/covariates.rds")

COVARS_LIST <- c("baseline_age", "age_sq", "FUyrs", "FU_n", 
                 "baseline_trait", "UKB_assmt_centre")

# Keep ids that pass QC and relevant covariates ----

full_dat <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    res <- blups[[p]][[sx]]
    # Subset to ids that passed QC
    res <- subset(res, res$eid %in% qcd_ids$IID)
    # Merge in covariates
    res <- merge(res, covars[[p]][, c("eid", "sex", 
                                      COVARS_LIST[!COVARS_LIST == "UKB_assmt_centre"])],
                 by = "eid")
    # Merge in UKB assessment centre
    res <- merge(res, qcd_ids[, c("IID", "UKB_assmt_centre")], 
                 by.x = "eid", by.y = "IID")
    return (res)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(full_dat) <- PHENOTYPES

# Adjust coefficients for covariates and RINT ----

adj_rint_slopes <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    df <- full_dat[[p]][[sx]]
    # For sex-combined strata, adjust for sex
    if (sx == "sex_comb") {
      COVARS_LIST <- c(COVARS_LIST, "sex")
    }
    # Formula for adjustment
    mod_formula <- 
      formula(paste0("lmm_slope ~ ", paste(COVARS_LIST, 
                                           collapse = " + ")))
    # Get residuals
    mod_resid <- lm(mod_formula, data = df)$residuals
    # RINT
    rinted_resid <- qnorm((rank(mod_resid) - 0.5) / sum(!is.na(mod_resid)))
    # Return formatted data
    res <- data.frame(FID = df$eid, 
                      IID = df$eid,
                      adj_slope = rinted_resid) 
    
    # Write results to table
    write.table(res,
                paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/traits_for_gwas/lmm_slopes_", 
                       p, "_", sx, ".txt"),
                sep = "\t", row.names = F, quote = F)
    return (res)
  })
})

