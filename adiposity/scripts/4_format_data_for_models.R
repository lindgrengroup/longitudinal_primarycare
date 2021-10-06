# Author: Samvida S. Venkatesh
# Date: 05/10/21

library(tidyverse)

# Read data ----

adipo <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/TRAINING_SET_adiposity.rds")
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/covariates.rds")

PHENOTYPES <- names(adipo)

# Combine covariates with adiposity data and split by sex ----

split_data <- lapply(PHENOTYPES, function (p) {
  full <- merge(adipo[[p]], covars[[p]], by = "eid")
  split_res <- split(full, full$sex)
  split_res[["sex_comb"]] <- full
  return (split_res)
})
names(dat) <- PHENOTYPES

saveRDS(split_data, 
        "/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/TRAINING_SET_adiposity_split_sex.rds")