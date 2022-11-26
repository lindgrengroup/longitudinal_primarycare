# Author: Samvida S. Venkatesh
# Date: 11/10/21

library(tidyverse)

# Read data ----

STRATA <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/strata_filenames.txt")$V1

trait_dat <- lapply(STRATA, function (strt) {
  print(paste0("Running strata: ", strt))
  get_files <- list.files("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/traits_for_GWAS/",
                          pattern = paste0("lmm_.*_", strt, ".txt"))
  
  res <- lapply(get_files, function (fname) {
    df <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/traits_for_GWAS/",
                      fname), header = T, stringsAsFactors = F)
    # Change colname (adj_trait) to trait name (slope, intercept, etc.)
    colnames(df)[3] <- gsub(paste0("_", strt, "*.txt"), "", fname)
    return (df)
  })
  
  # Bind columns (merge) for all traits, i.e. intercept, slope, etc.
  res <- res %>% reduce(inner_join, by = c("FID", "IID"))

  # Merge in covariate information from ids that passed unrelated-sample-QC
  cov_qc <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/sample_qc/", 
                              strt, "_UNRELATED_ids_passed_qc.txt"),
                       sep = "\t", header = T, stringsAsFactors = F)
  res <- res %>% inner_join(cov_qc, by = c("FID", "IID"))
  res <- res[complete.cases(res), ]
  
  # Write data
  write.table(res,
              paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/traits_for_GWAS/lmm_traits_",
                     strt, ".txt"),
              sep = "\t", quote = F, row.names = F)
  
  return (res)
})
