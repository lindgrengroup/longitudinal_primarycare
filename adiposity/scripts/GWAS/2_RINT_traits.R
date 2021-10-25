# Author: Samvida S. Venkatesh
# Date: 11/10/21

library(tidyverse)

# Read data ----

PHENOTYPES <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/pheno_names.txt",
                         sep = "\t", header = F, stringsAsFactors = F)$V1

# Random effect coefficients
rand_effs <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/random_effect_terms_",
                 p, ".rds"))
})
names(rand_effs) <- PHENOTYPES
SEX_STRATA <- names(rand_effs[[1]])

# IDs that passed sample QC for GWAS 
qcd_ids <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/ids_passed_qc_211015.txt", 
                      sep = "\t", header = T)

# Covariates file (general and trait-specific) 
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/covariates.rds")

COVARS_LIST <- c("baseline_trait", "UKB_assmt_centre")

# Keep ids that pass QC and relevant covariates ----

full_dat <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    res <- rand_effs[[p]][[sx]]
    # Subset to ids that passed QC
    res <- subset(res, res$eid %in% qcd_ids$IID)
    # Remove (intercept) column as that's not used for GWAS
    res <- res[, which(colnames(res) != "(Intercept)")]
    # Merge in covariates
    res <- merge(res, covars[[p]][, c("eid", "sex", "baseline_trait")],
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

adj_and_rint <- function (tm_name, sx, dat) {
  # For sex-combined strata, adjust for sex
  if (sx == "sex_comb") {
    COVARS_LIST <- c(COVARS_LIST, "sex")
  }
  # Formula for adjustment
  mod_formula <- 
    formula(paste0(tm_name, " ~ ", paste(COVARS_LIST, 
                                         collapse = " + ")))
  # Get residuals
  mod_resid <- lm(mod_formula, data = dat)$residuals
  # RINT
  rinted_resid <- qnorm((rank(mod_resid) - 0.5) / sum(!is.na(mod_resid)))
  # Return formatted data
  res <- data.frame(eid = dat$eid, adj_tm = rinted_resid)
  colnames(res)[2] <- paste0("adj_", tm_name)
  return (res)
}

adj_rint_slopes <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    df <- full_dat[[p]][[sx]]
    tms_to_rint <- colnames(df)[!colnames(df) %in% 
                                  c("eid", "IID", "sex", COVARS_LIST)]
    rinted_cols <- lapply(tms_to_rint, function (tm) {
      return (adj_and_rint(tm, sx, df))
    })
    rinted_cols <- rinted_cols %>% reduce(inner_join, by = "eid") %>%
      mutate(FID = eid, IID = eid)
    rinted_cols <- rinted_cols[, c("FID", "IID",
                                   colnames(rinted_cols)[grep("^adj", 
                                                              colnames(rinted_cols))])]
    # Write results to table
    write.table(rinted_cols,
                paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/traits_for_gwas/", 
                       p, "_", sx, ".txt"),
                sep = "\t", row.names = F, quote = F)
    return (rinted_cols)
  })
})

