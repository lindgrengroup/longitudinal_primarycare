# Author: Samvida S. Venkatesh
# Date: 11/10/21

library(tidyverse)

# Read data ----

# Name of term to RINT 
TERM_TO_RINT <- "ns.*1"
result_prefix <- "/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/traits_for_gwas/cspline_ns1"

# BLUP files
PHENOTYPES <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/pheno_names.txt")$V1
blups <- lapply(PHENOTYPES, function (p) {
  res <- readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/without_pcs_cubic_time_spline_blups_",
                        p, ".rds"))
  res <- lapply(res, function (df) {
    eids <- rownames(df)
    model_terms <- df[, grep(TERM_TO_RINT, colnames(df))]
    to_return <- data.frame(eid = eids,
                            model_term = model_terms)
    return (to_return)
  })
  return (res)
})
names(blups) <- PHENOTYPES

SEX_STRATA <- names(blups[[1]])

# IDs that passed sample QC
ids_passed_qc <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SEX_STRATA, function (sx) {
    read.table(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/sample_qc/", 
                      p, "_", sx, "_ids_passed_qc.txt"),
               sep = "\t", header = T)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(ids_passed_qc) <- PHENOTYPES

# Covariates file (general and trait-specific) 
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/covariates.rds")[PHENOTYPES]

COVARS_LIST <- c("baseline_age", "age_sq", "FUyrs", "FU_n", 
                 "UKB_assmt_centre")

# Extract BLUP terms and relevant covariates ----

full_dat <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    # Filter to samples that passed QC (previously calculated)
    res <- blups[[p]][[sx]] %>% 
      filter(eid %in% ids_passed_qc[[p]][[sx]]$IID)
    # Merge in covariates
    res <- merge(res, covars[[p]][, c("eid", "sex", 
                                      COVARS_LIST[!COVARS_LIST == "UKB_assmt_centre"])],
                 by = "eid")
    # Merge in UKB assessment centre
    res <- merge(res, 
                 ids_passed_qc[[p]][[sx]][, c("IID", "UKB_assmt_centre")], 
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
      formula(paste0("model_term ~ ", paste(COVARS_LIST, 
                                           collapse = " + ")))
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
                paste0(result_prefix, "_", p, "_", sx, ".txt"),
                sep = "\t", row.names = F, quote = F)
    return (res)
  })
  return ()
})

