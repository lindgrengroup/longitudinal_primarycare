# Author: Samvida S. Venkatesh
# Date: 07/02/22

library(tidyverse)
library(nnet)

# Read data ----

result_prefix <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/traits_for_GWAS/"

# Clustering results 
PHENOTYPES <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/pheno_names.txt")$V1
SEX_STRATA <- c("F", "M", "sex_comb")

cluster_results <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/assigned_clusters_",
                             p, "_", sx, ".txt"),
                      sep = "\t", header = T, stringsAsFactors = F)
    res$eid <- as.character(res$eid)
    res <- res[, c("eid", "clust")]
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

# Covariates file (general and trait-specific) 
general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[PHENOTYPES]
covars <- lapply(covars, function (df) {
  df$eid <- as.character(df$eid)
  return (df)
})

PCs <- paste0("PC", c(1:21))
COVARS_LIST <- c("UKB_assmt_centre", "genotyping.array", "year_of_birth",
                 "baseline_age", "age_sq", "baseline_trait", "sex")

# Add covariates ----

full_dat <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    # Filter to samples that passed QC (previously calculated)
    res <- cluster_results[[p]][[sx]] %>% 
      filter(eid %in% ids_passed_qc[[p]][[sx]]$IID) %>%
      mutate(clust = factor(as.character(clust)))
    
    # Merge in covariates
    res <- left_join(res, general_covars %>% select(any_of(c("eid", PCs, COVARS_LIST))), 
                     by = "eid")
    res <- left_join(res, covars[[p]] %>% select(any_of(c("eid", COVARS_LIST))),
                     by = "eid")
    res$genotyping.array <- ids_passed_qc[[p]][[sx]]$genotyping.array[match(res$eid,
                                                                            ids_passed_qc[[p]][[sx]]$IID)]
    return (res)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(full_dat) <- PHENOTYPES

# Write variance explained by covariates ----

getLogLik <- function (df, cov_list) {
  var_form <- formula(paste0("clust ~ ", paste0(cov_list, collapse = " + ")))
  sub_df <- df[, c("eid", "clust", cov_list)]
  sub_df <- sub_df[complete.cases(sub_df), ]
  
  mod_ll <- logLik(multinom(var_form, data = sub_df))
  return (mod_ll)
} 

var_sets_to_test <- list(PCs = PCs, 
                         sex = "sex",
                         age = "baseline_age", 
                         age_age_sq = c("baseline_age", "age_sq"),
                         birth_year = "year_of_birth",
                         age_birth_year = c("baseline_age", "age_sq", "year_of_birth"),
                         all_covs = c(PCs, "sex", "baseline_age", "age_sq",
                                      "year_of_birth", "UKB_assmt_centre",
                                      "genotyping.array"))

# var_sets_to_test <- list(PCs = PCs, 
#                          sex = "sex",
#                          age = "baseline_age", 
#                          age_age_sq = c("baseline_age", "age_sq"),
#                          birth_year = "year_of_birth",
#                          age_birth_year = c("baseline_age", "age_sq", "year_of_birth"),
#                          baseline_trait = "baseline_trait",
#                          all_covs = c(PCs, "sex", "baseline_age", "age_sq",
#                                       "year_of_birth", "UKB_assmt_centre",
#                                       "genotyping.array", "baseline_trait"))

var_set_names <- names(var_sets_to_test)

var_expl_table <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    df <- full_dat[[p]][[sx]]
    null_ll <- logLik(multinom(clust ~ 1, data = df))
    
    var_lists <- lapply(var_set_names, function (vn) {
      pseudo_r2 <- 1 - (getLogLik(df, var_sets_to_test[[vn]]) / null_ll)
      res <- data.frame(pheno = p, sex_strata = sx,
                        pseudo_r2 = pseudo_r2,
                        vars_tested = vn)
      return (res)
    })
    var_lists <- bind_rows(var_lists)
    write.table(var_lists,
                paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/results/var_explained_", 
                       p, "_", sx, ".txt"),
                sep = "\t", row.names = F, quote = F)
    return ()
  })
  return ()
})
return ()

# Write wide form data ----

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    
    # Wrangle data to wide form (k1, k2, etc. for cluster belonging)
    res <- full_dat[[p]][[sx]] %>% mutate(clust = paste0("k", clust)) %>%
      pivot_wider(id_cols = -clust, 
                  names_from = clust, 
                  values_from = clust, 
                  values_fn = function (x) 1, 
                  values_fill = 0)
    
    # Create group cols
    res <- res %>% 
      mutate(k1_k2 = ifelse(k1 + k2 > 0, 1, 0),
             k1_k2_k3 = ifelse(k1 + k2 + k3 > 0, 1, 0))
    
    # Write results to table
    write.table(res,
                paste0(result_prefix, "cluster_membership_", p, "_", sx, ".txt"),
                sep = " ", row.names = F, quote = F)
    
    return ()
  })
  return ()
})