# Author: Samvida S. Venkatesh
# Date: 11/10/21

library(tidyverse)

# Read data ----

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

ADJ_COVARS <- c("baseline_age", "age_sq", "year_of_birth", "FU_n", "FUyrs",
                "UKB_assmt_centre") # add sex for sc analyses

blups <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SEX_STRATA, function (sx) {
    df <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/lmm_models/",
                            p, "_", sx, "_all_blups.txt"), 
                     sep = "\t", header = T, stringsAsFactors = F)
    df$eid <- as.character(df$eid)
    return (df)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(blups) <- PHENOTYPES

clustprobs <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SEX_STRATA, function (sx) {
    df <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/highdim_splines_clustering/", 
                            p, "_", sx, "/soft_clustering_probs_", p, "_", sx, ".txt"), 
                     sep = "\t", header = T, stringsAsFactors = F)
    df <- df %>% mutate(eid = as.character(eid),
                        k1_k2 = k1 + k2,
                        k1_k2_k3 = k1 + k2 + k3)
    return (df)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(clustprobs) <- PHENOTYPES

# Raw data (original) to calculate covariates
raw_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/data/main_data_adipo_change.rds")[PHENOTYPES]

# Covariates
general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

# IDs that passed sample QC
ids_passed_qc <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SEX_STRATA, function (sx) {
    df <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/replication_genotypes/sample_qc/", 
                            p, "_", sx, "_ids_passed_qc.txt"),
                     sep = "\t", header = T)
    df$IID <- as.character(df$IID)
    df$FID <- as.character(df$FID)
    return (df)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(ids_passed_qc) <- PHENOTYPES

# Get covariates for each phenotype and strata ----

covars <- lapply(PHENOTYPES, function (p) {
  res <- raw_dat[[p]] %>%
    group_by(eid) %>% 
    arrange(age_event, .by_group = T) %>%
    summarise(baseline_trait = first(value),
              baseline_age = first(age_event),
              age_sq = baseline_age^2,
              FU_n = n(),
              FUyrs = last(age_event) - first(age_event))
  res <- left_join(res, 
                   general_covars[, c("eid", "sex", "year_of_birth",
                                      "UKB_assmt_centre")],
                   by = "eid")
  return (res)
})
names(covars) <- PHENOTYPES

# Functions to transform the probability traits ----

squeezeProbs <- function (x, nsamples = 100) {
  squeezed_x <- (x*(nsamples - 1) + 0.5)/nsamples
  return (squeezed_x)
}

getLogit <- function (x) {
  return (log(x/(1-x)))
}

# Create phenotype files ----

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    print(paste0(p, "_", sx))
    lmm_df <- blups[[p]][[sx]]
    # Only retain the IDs that passed genotyping sample QC
    lmm_df <- lmm_df %>% filter(eid %in% ids_passed_qc[[p]][[sx]]$IID)
    # Add covariates
    lmm_df <- left_join(lmm_df, covars[[p]], by = "eid") 
    
    # Linearly adjust slope for covariates
    to_adj <- ADJ_COVARS
    if (sx == "sex_comb") to_adj <- c(ADJ_COVARS, "sex")
    model_formula <- formula(paste0("t ~ X.Intercept. + ", paste0(ADJ_COVARS, collapse = " + ")))
    
    modelled_slope <- lm(model_formula, lmm_df)
    # RINT
    trait_to_rint <- resid(modelled_slope)
    rinted_trait <- qnorm((rank(trait_to_rint) - 0.5) / sum(!is.na(trait_to_rint)))
    
    # Lmm phenotypes
    pheno_dat <- data.frame(eid = lmm_df$eid,
                            b1 = rinted_trait)
    
    # Soft clustering phenotypes
    logit_probs <- as.data.frame(apply(clustprobs[[p]][[sx]][, c("k1", "k1_k2", "k1_k2_k3")], 2, 
                         FUN = function (x) {
                           getLogit(squeezeProbs(x, nsamples = 100))
                         }))
    logit_probs$eid <- clustprobs[[p]][[sx]]$eid
    pheno_dat <- left_join(pheno_dat, logit_probs, 
                           by = "eid")
    
    # Gather covariates for PLINK file, present in the sample QC file
    to_write <- inner_join(pheno_dat, covars[[p]])
    # Change sex coding to 1/0 for PLINK: 1 = male, 0 = female
    to_write$sex <- ifelse(to_write$sex == "M", 1, 0)
    id_info <- ids_passed_qc[[p]][[sx]] %>% select(-any_of(c("UKB_assmt_centre")))
    to_write <- inner_join(id_info, 
                           to_write, 
                           by = c("IID" = "eid"))
    # Change colname to start with "#" for PLINK
    colnames(to_write)[1] <- "#FID"
    write.table(to_write,
                paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/replication_genotypes/traits_for_GWAS/",
                       p, "_", sx, ".txt"),
                sep = "\t", quote = F, row.names = F)
  })
})
