# Author: Samvida S. Venkatesh
# Date: 07/02/22

library(tidyverse)

# Read data ----

result_prefix <- "/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/traits_for_gwas/"

# Clustering results 
PHENOTYPES <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/pheno_names.txt")$V1
SEX_STRATA <- c("F", "M", "sex_comb")

cluster_results <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    res <- read.table(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/age_adj_clusters_",
                             p, "_", sx, ".txt"),
                      sep = "\t", header = T, stringsAsFactors = F)
    return (res)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(cluster_results) <- PHENOTYPES

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

PCs <- paste0("PC", c(1:21))
COVARS_LIST <- c("UKB_assmt_centre", "genotyping.array")

# Wrangle cluster data into wide form and get relevant covariates ----

full_dat <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    # Filter to samples that passed QC (previously calculated)
    res <- cluster_results[[p]][[sx]] %>% 
      filter(eid %in% ids_passed_qc[[p]][[sx]]$IID)
    
    # Wrangle data to wide form (k1, k2, etc. for cluster belonging)
    res <- res %>% pivot_wider(id_cols = eid, 
                               names_from = k, names_prefix = "k",
                               values_from = k, values_fn = function (x) 1, 
                               values_fill = 0)
    
    # Merge in sex 
    res <- merge(res, covars[[p]][, c("eid", "sex")],
                 by = "eid")
    # Merge in UKB assessment centre and genotyping array
    res <- merge(res, 
                 ids_passed_qc[[p]][[sx]][, c("IID", COVARS_LIST, PCs)], 
                 by.x = "eid", by.y = "IID")
    
    # Write results to table
    write.table(res,
                paste0(result_prefix, p, "_", sx, "/cluster_membership_",
                p, "_", sx, ".txt"),
                sep = " ", row.names = F, quote = F)
    
    return ()
  })
 return ()
})