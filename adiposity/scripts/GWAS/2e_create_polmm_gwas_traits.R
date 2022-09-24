# Author: Samvida S. Venkatesh
# Date: 15/09/2022

library(tidyverse)

# Read data ----

result_prefix <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/traits_for_GWAS/"

# Clustering results 
PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

cluster_results <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/",
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
  return (res)
})

# Get cluster assignment based on max probability ----

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    # Only keep samples that pass QC
    res <- cluster_results[[p]][[sx]]
    res <- res %>% filter(eid %in% ids_passed_qc[[p]][[sx]]$IID)
    res$genotyping.array <- ids_passed_qc[[p]][[sx]]$genotyping.array[match(res$eid,
                                                                            ids_passed_qc[[p]][[sx]]$IID)]
    
    # Get cluster assignment
    for_assn <- res %>% select(-all_of(c("eid", "genotyping.array")))
    clust_assigned <- apply(for_assn, 1, 
                            function (x) colnames(for_assn)[which.max(x)])
    
    # Return formatted data
    to_write <- data.frame(eid = res$eid, 
                           clust = clust_assigned,
                           genotyping.array = res$genotyping.array) 
    # Add in covariates
    to_write <- left_join(to_write, 
                    covars[[p]], by = "eid")
    
    write.table(to_write, 
                paste0(result_prefix, p, "_", sx, "_assigned_clusts.txt"),
                sep = "\t", row.names = F, quote = F)
  })
})
