# Author: Samvida S. Venkatesh
# Date: 05/10/21

library(tidyverse)

set.seed(051021)

# Read adiposity and covariates data ----

adiposity <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/QCd_adiposity.rds")
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/covariates.rds")

PHENOTYPES <- names(adiposity)

# Split IDs into test and training sets ----

TEST_FRACTION <- 0.2

test_ids <- lapply(covars, function (df) {
  all_ids <- unique(df$eid)
  test_ids <- sample(all_ids, size = TEST_FRACTION*length(all_ids), replace = F)
  return (test_ids)
})
names(test_ids) <- PHENOTYPES

test_adiposity <- lapply(PHENOTYPES, function (p) {
  res <- adiposity[[p]]
  res <- subset(res, res$eid %in% test_ids[[p]])
  # only retain columns of interest
  res <- res[, c("eid", "data_provider", "event_dt", "age_event",
                 "value", "biomarker")]
  return (res)
})
names(test_adiposity) <- PHENOTYPES

train_adiposity <- lapply(PHENOTYPES, function (p) {
  res <- adiposity[[p]]
  res <- subset(res, !res$eid %in% test_ids[[p]])
  res <- res[, c("eid", "data_provider", "event_dt", "age_event",
                 "value", "biomarker")]
  return (res)
})
names(train_adiposity) <- PHENOTYPES

# Sanity check that there is no overlap between test and train sets ----

lapply(PHENOTYPES, function (p) {
  overlap <- inner_join(test_adiposity[[p]], train_adiposity[[p]])
  if (nrow(overlap) == 0) {
    cat(paste0("** PHENOTYPE **", p, "\n", 
               "No overlap between test and train sets", "\n"))
  } else { cat(paste0("** PHENOTYPE **", p, "\n", 
                      "*** STOP ERROR: OVERLAP BETWEEN TEST AND TRAIN SETS ***",
                      "\n")) }
  return ()
})

# Save test and train sets ----

saveRDS(test_adiposity, 
        "/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/TEST_SET_adiposity.rds")
saveRDS(train_adiposity, 
        "/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/TRAINING_SET_adiposity.rds")



