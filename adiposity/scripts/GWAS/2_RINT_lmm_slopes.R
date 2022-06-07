# Author: Samvida S. Venkatesh
# Date: 11/10/21

library(tidyverse)

# Read data ----

# BLUP files
PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

blups <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SEX_STRATA, function (sx) {
    df <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/lmm_models/",
                            p, "_", sx, "_blups_full_model.txt"),
                     sep = "\t", header = T, stringsAsFactors = F)
    return (df)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(blups) <- PHENOTYPES

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

# Create RINTed slope adjusted for intercept ----

rinted_slope_adj_int <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    df <- blups[[p]][[sx]]
    # Only retain the IDs that passed genotyping sample QC
    df <- df %>% filter(eid %in% ids_passed_qc[[p]][[sx]]$IID)
    
    # Linearly adjust slope for intercept
    modelled_slope <- lm(t ~ X.Intercept., df)
    
    # RINT
    trait_to_rint <- resid(modelled_slope)
    rinted_trait <- qnorm((rank(trait_to_rint) - 0.5) / sum(!is.na(trait_to_rint)))
    # Return formatted data
    res <- data.frame(FID = df$eid, 
                      IID = df$eid,
                      rinted_trait = rinted_trait) 
    
    # Write results to table
    write.table(res,
                paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/traits_for_GWAS/lmm_slopes_adj_int_", 
                       p, "_", sx, ".txt"),
                sep = "\t", row.names = F, quote = F)
    return ()
  })
  return ()
})

# Create RINTed slope within each half of intercept ----

rinted_slope_halves <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    df <- blups[[p]][[sx]]
    # Only retain the IDs that passed genotyping sample QC
    df <- df %>% filter(eid %in% ids_passed_qc[[p]][[sx]]$IID)
    
    # Create top and bottom halves of data by intercept
    med_intercept <- median(df$X.Intercept.)
    top_int <- df %>% filter(X.Intercept. >= med_intercept)
    bottom_int <- df %>% filter(X.Intercept. < med_intercept)
    
    # RINT within each half
    
    # Top half
    trait_to_rint <- top_int$t
    rinted_trait <- qnorm((rank(trait_to_rint) - 0.5) / sum(!is.na(trait_to_rint)))
    # Return formatted data
    res <- data.frame(FID = df$eid, 
                      IID = df$eid,
                      rinted_trait = rinted_trait)
    # Write results to table
    write.table(res,
                paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/traits_for_GWAS/lmm_slopes_top_int_", 
                       p, "_", sx, ".txt"),
                sep = "\t", row.names = F, quote = F)
    
    # Bottom half
    trait_to_rint <- bottom_int$t
    rinted_trait <- qnorm((rank(trait_to_rint) - 0.5) / sum(!is.na(trait_to_rint)))
    # Return formatted data
    res <- data.frame(FID = df$eid, 
                      IID = df$eid,
                      rinted_trait = rinted_trait)
    # Write results to table
    write.table(res,
                paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/traits_for_GWAS/lmm_slopes_bottom_int_", 
                       p, "_", sx, ".txt"),
                sep = "\t", row.names = F, quote = F)
  
    return ()
  })
  return ()
})
