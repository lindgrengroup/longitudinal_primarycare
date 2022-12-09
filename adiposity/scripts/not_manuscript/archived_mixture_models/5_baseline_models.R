# Author: Samvida S. Venkatesh
# Date: 19/05/2021

library(tidyverse)
theme_set(theme_bw())
library(flexmix)

set.seed(190521)

# Read files ----

adiposity <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/visually_QCd_adiposity.rds")
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/model_covariates.rds")
PCs <- paste0("PC", 1:21)

# Read arguments and prepare subsets for analysis ----

args <- commandArgs(trailingOnly = T)

PHENO <- args[1]
STRATUM <- args[2]
STRATUM_ANCESTRY <- sub("_.*", "", STRATUM)
STRATUM_SEX <- sub(".*_", "", STRATUM)

sub_covars <- covars[[PHENO]]
adipo <- adiposity[[PHENO]]

# Subset adiposity data relevant to stratum
# One last check to remove any individuals without multiple measures

si_eids <- sub_covars$eid[sub_covars$ancestry == STRATUM_ANCESTRY &
                            sub_covars$FU_n > 1]
if(STRATUM_SEX != "sexcomb") {
  si_eids <- sub_covars$eid[sub_covars$sex == STRATUM_SEX &
                              sub_covars$ancestry == STRATUM_ANCESTRY &
                              sub_covars$FU_n > 1]
}
dat <- subset(adipo, adipo$eid %in% si_eids)

# Add covariates to the data to model
dat <- merge(dat, sub_covars, by = "eid", all.x = T)

# Adjust value for covariates
adj_formula <- ifelse(STRATUM_SEX == "sexcomb", 
                      paste0("value ~ baseline_age + age_sq + baseline_trait + FUyrs + FU_n + sex + ",
                             paste(PCs, collapse = " + ")),
                      paste0("value ~ baseline_age + age_sq + baseline_trait + FUyrs + FU_n + ",
                             paste(PCs, collapse = " + "))) 
adj_value_model <- lm(formula(adj_formula), data = dat)
dat$adj_value <- adj_value_model$residuals

saveRDS(dat, paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/adj_traits_and_covars/", 
                    PHENO, "_", STRATUM, ".rds"))

# Baseline models (1 group) ----

# Linear
mod_LIN <- FLXMRlmm(adj_value ~ 0 + age_event, 
                    random = ~ 1 + age_event,
                    varFix = c(Random = F, Residual = F))
# Smoothing spline 
mod_SPL <- FLXMRlmm(adj_value ~ 0 + age_event, 
                    random = ~ 1 + age_event, 
                    lm.fit = "smooth.spline",
                    varFix = c(Random = F, Residual = F))

# Flexmix for baseline models

runFlex <- function (mod) {
  out <- tryCatch(
    expr = {
      # Run model
      res <- suppressWarnings(flexmix(.~.|eid, k = 1, model = mod, 
              data = dat, 
              control = list(iter.max = 1000, minprior = 0)))
      message("Running baseline model...")
      res
    }, 
    error = function (cond) {
      message("Model fails to run, original error message:")
      message(cond)
      # Return value in case of error
      return (NA)
    }
  )
  return (out)
}

flex_bl_LIN <- runFlex(mod_LIN)
flex_bl_SPL <- runFlex(mod_SPL)

if (!is.na(flex_bl_LIN)) { 
  print(summary(flex_bl_LIN)) 
  print ("\n") }

if (!is.na(flex_bl_SPL)) { 
  print(summary(flex_bl_SPL)) 
  print ("\n") }

# Report metrics
sink(paste0("log_files/mixture_models/baseline_models_", PHENO, ".txt"), 
     append = T)
cat(paste0("Strata: ", 
           STRATUM_ANCESTRY, " ", STRATUM_SEX, "\n"))
if (!is.na(flex_bl_LIN)) { 
  print(summary(flex_bl_LIN)) 
  print ("\n") }
if (!is.na(flex_bl_SPL)) { 
  print(summary(flex_bl_SPL)) 
  print ("\n") }
sink()
