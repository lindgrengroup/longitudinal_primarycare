# Author: Samvida S. Venkatesh
# Date: 11/10/21

library(tidyverse)

# Read data ----

all_files <- list.files("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/traits_for_gwas/",
                        pattern = "*.txt")
TRAIT_NAMES <- gsub("*.txt", "", all_files)

trait_dat <- lapply(TRAIT_NAMES, function (trait) {
  read.table(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/traits_for_gwas/", 
                    trait, ".txt"), header = T, stringsAsFactors = F)
})
names(trait_dat) <- TRAIT_NAMES

# Get sample order from UKB sample file
sample_order <- read.table("/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample",
                           sep = " ", header = T, stringsAsFactors = F)$ID_1

# Wrangle data into correct sample order -----

ordered_dat <- lapply(TRAIT_NAMES, function (trait) {
  df <- trait_dat[[trait]]
  res <- data.frame(IID = sample_order)
  # Order trait by sample file and replace missing values with NA
  res <- left_join(res, df[, c("IID", "adj_trait")],
                   by = "IID") %>%
    replace_na(list(adj_trait = -999))
  # Only return a single numeric vector with the trait values
  return (as.numeric(res$adj_trait))
})
ordered_dat <- bind_cols(ordered_dat)
colnames(ordered_dat) <- TRAIT_NAMES

write.table(ordered_dat,
            "/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/traits_for_gwas/all_traits_for_bgenie.txt",
            sep = " ", quote = F, row.names = F)
