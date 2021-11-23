# Author: Samvida S. Venkatesh
# Date: 11/10/21

library(tidyverse)

# Read data ----

# Name of term to RINT 
TERM_TO_RINT <- "t"
result_prefix <- "/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/traits_for_gwas/lmm_slopes_no_adj"

# BLUP files
PHENOTYPES <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/pheno_names.txt")$V1
blups <- lapply(PHENOTYPES, function (p) {
  res <- readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/lmm_blups_",
                        p, ".rds"))
  res <- lapply(res, function (df) {
    b0s <- df[, "(Intercept)"]
    b0_quartile <- cut(b0s, 
                       breaks = c(quantile(b0s, probs = seq(0, 1, by = 1/4))),
                       labels = c(1:4),
                       include.lowest = T)
    
    df$eid <- rownames(df)
    df <- df[, c("eid", TERM_TO_RINT)]
    colnames(df) <- c("eid", "model_term")
    df$b0_quartile <- b0_quartile
    return (df)
  })
  return (res)
})
names(blups) <- PHENOTYPES

SEX_STRATA <- names(blups[[1]])
NPCs <- 21

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

COVARS_LIST <- c("baseline_age", "age_sq", "UKB_assmt_centre")

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

# Within each quartile, adjust coefficients for covariates and RINT ----

adj_rint_slopes <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    # For sex-combined strata, adjust for sex
    if (sx == "sex_comb") {
      COVARS_LIST <- c(COVARS_LIST, "sex")
    }
    # Formula for adjustment
    mod_formula <- 
      formula(paste0("model_term ~ ", paste(COVARS_LIST, 
                                            collapse = " + ")))
    
    # Do this separately in each quartile of b0 
    res <- lapply(c(1:4), function (q) {
      dat_to_model <- full_dat[[p]][[sx]] %>% filter(b0_quartile == q)
      # Get residuals
      mod_resid <- lm(mod_formula, data = dat_to_model)$residuals
      # RINT
      rinted_resid <- qnorm((rank(mod_resid) - 0.5) / sum(!is.na(mod_resid)))
      # Return formatted data
      to_write <- data.frame(FID = dat_to_model$eid, 
                        IID = dat_to_model$eid,
                        adj_trait = rinted_resid) 
      
      # Write results to table
      write.table(to_write,
                  paste0(result_prefix, "_b0_quartile_", q, "_", 
                         p, "_", sx, ".txt"),
                  sep = "\t", row.names = F, quote = F)
    })
    return (res)
  })
  return ()
})

