# Author: Samvida S. Venkatesh
# Date: 21/10/21

library(lme4)
library(tidyverse)

# Read data ----

models <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/lmm_all_strata.rds")
PHENOTYPES <- names(models)
SEX_STRATA <- names(models[[1]])

covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/covariates.rds")

newdat <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/TEST_SET_adiposity.rds")

log_file <- "/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/log_files/lmm_prediction_accuracy.txt"

# Add covariates and time from baseline to new data ----

newdat <- lapply(PHENOTYPES, function (p) {
  res <- merge(newdat[[p]], covars[[p]], by = "eid")
  return (res %>% mutate(t = age_event - baseline_age))
})
names(newdat) <- PHENOTYPES

# Predict from models ----

pred_values <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    
    if (sx == "sex_comb") {
      newdf <- newdat[[p]]
    } else {
      newdf <- newdat[[p]] %>% filter(sex == sx)
    } 
    
    pred_res <- as.data.frame(predict(models[[p]][[sx]],
                        newdata = newdf, allow.new.levels = T))
    colnames(pred_res) <- "predicted_value"
    full_res <- bind_cols(newdf, pred_res)
    return (full_res)
  })
  names(per_sex) <- SEX_STRATA
  return (per_sex)
})
names(pred_values) <- PHENOTYPES

saveRDS(pred_values, 
        "/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/TEST_SET_lmm_predictions.rds")

# Print MAE and MSE for all strata ----

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    errors <- mutate(pred_values[[p]][[sx]],
                     error = predicted_value - value,
                     sqerr = error^2,
                     abserr = abs(error))
    sink(log_file, append = T)
    cat(paste0("** Phenotype and Strata: ", p, ", ", sx, "\n",
               "\t", "MSE = ", mean(errors$sqerr, na.rm = T), "\n",
               "\t", "MAE = ", mean(errors$abserr, na.rm = T), "\n"))
    sink()
  })
})
