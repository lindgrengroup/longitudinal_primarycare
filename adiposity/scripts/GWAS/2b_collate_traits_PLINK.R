# Author: Samvida S. Venkatesh
# Date: 11/10/21

library(tidyverse)

# Read data ----

STRATA <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/strata_filenames.txt")$V1

trait_dat <- lapply(STRATA, function (strt) {
  get_files <- list.files("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/traits_for_gwas/",
                          pattern = paste0(strt, "*.txt"))
  
  res <- lapply(get_files, function (fname) {
    df <- read.table(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/traits_for_gwas/",
                      fname), header = T, stringsAsFactors = F)
    # Change colname (adj_trait) to trait name (slope, intercept, etc.)
    colnames(df)[3] <- gsub(paste0("_", strt, "*.txt"), "", fname)
    return (df)
  })
  names(res) <- gsub(paste0("_", strt, "*.txt"), "", get_files)
  
  # Bind columns (merge) for all traits, i.e. intercept, slope, etc.
  res <- res %>% reduce(full_join, by = c("FID", "IID"))
  # Change colname to start with "#" for PLINK
  colnames(res)[1] <- "#FID"
  
  # Write data
  write.table(res,
              paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/traits_for_gwas/all_models_",
                     strt, ".txt"),
              sep = "\t", quote = F, row.names = F)
  
  return (res)
})
names(trait_dat) <- STRATA
