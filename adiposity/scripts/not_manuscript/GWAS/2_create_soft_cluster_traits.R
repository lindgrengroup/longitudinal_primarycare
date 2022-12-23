# Author: Samvida S. Venkatesh
# Date: 24/08/22

library(tidyverse)

# Read data ----

result_prefix <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/GWAS/traits_for_GWAS/"

# Clustering results 
PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

cluster_results <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/clustering/",
                             p, "_", sx, "/soft_clustering_probs_", p, "_", sx, ".txt"),
                      sep = "\t", header = T, stringsAsFactors = F)
    res$eid <- as.character(res$eid)
    return (res)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(cluster_results) <- PHENOTYPES

# IDs that passed sample QC
ids_passed_qc <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/sample_qc/", 
                             p, "_", sx, "_ids_passed_qc.txt"),
                      sep = "\t", header = T)
    res$IID <- as.character(res$IID)
    return (res)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(ids_passed_qc) <- PHENOTYPES

# Covariates
covars <- lapply(PHENOTYPES, function (p) {
  readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[[p]]
})
names(covars) <- PHENOTYPES

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

covars <- lapply(covars, function (df) {
  res <- merge(df, general_covars, by = "eid")
  res$FID <- as.character(res$eid)
  res$IID <- as.character(res$eid)
  return (res)
})

# Already contains genotyping array, UKB asst centre, and PCs
ADD_COVARS_LIST <- c("sex", "baseline_trait", "baseline_age", "age_sq", "FUyrs", "FU_n")

# Logit-transform probabilities to create phenotype for linear GWAS ----

# Based on https://psycnet.apa.org/record/2006-03820-004?doi=1
# Smithson & Verkuilen 2006
squeezeProbs <- function (x, nsamples = 100) {
  squeezed_x <- (x*(nsamples - 1) + 0.5)/nsamples
  return (squeezed_x)
}

getLogit <- function (x) {
  return (log(x/(1-x)))
}

## CHECK THIS IF RE-RUNNING
# SPECIFY NBOOTSTRAP SAMPLES
# SPECIFY NCLUSTERS

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    # Only keep samples that pass QC
    res <- cluster_results[[p]][[sx]]
    res <- res %>% filter(eid %in% ids_passed_qc[[p]][[sx]]$IID)
    
    # Create GWAS traits
    res <- res %>% 
      mutate(k1_k2 = k1 + k2,
             k1_k2_k3 = k1 + k2 + k3)
    
    logit_probs <- apply(res[, c("k1", "k1_k2", "k1_k2_k3")], 2, 
                         FUN = function (x) {
      getLogit(squeezeProbs(x, nsamples = 100))
    })
    
    # Return formatted data
    to_write <- data.frame(FID = res$eid, 
                           IID = res$eid,
                           k1 = logit_probs[, "k1"], 
                           k1_k2 = logit_probs[, "k1_k2"],
                           k1_k2_k3 = logit_probs[, "k1_k2_k3"]) 
    write.table(to_write, 
                paste0(result_prefix, p, "_", sx, "_soft_clust_probs.txt"),
                sep = "\t", row.names = F, quote = F)
  })
})

# Write covariates file with relevant covariates for GWAS ----

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    df <- ids_passed_qc[[p]][[sx]]
    df$FID <- as.character(df$FID)
    df$IID <- as.character(df$IID)
    df <- left_join(df, 
                    covars[[p]][, c("IID", "FID", ADD_COVARS_LIST)], by = c("IID", "FID"))
    write.table(df, 
                paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/GWAS/sample_qc/", 
                       p, "_", sx, "_covariates.txt"),
                sep = "\t", row.names = F, quote = F)
  })
})

