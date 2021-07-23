# Author: Samvida S. Venkatesh
# Date: 15/06/21

library(lme4)
library(splines)
library(merTools)
library(lmerTest)
library(tidyverse)
theme_set(theme_bw())

set.seed(150621)

# Parse in phenotype argument
args <- commandArgs(trailingOnly = T)
PHENO <- args[1]

SEX_STRATA <- c("F", "M", "sex_comb")

# Read files ----

dat <- 
  readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/adj_traits/adj_",
                 PHENO, ".rds"))

# Log file
log_file <- paste0("log_files/traj_modelling_log_", PHENO, ".txt")

# Plot basis functions for non-linearity ----

# Scale the age term to avoid multicollinearity
dat <- lapply(dat, function (df) {
  df$scaled_age <- scale(df$age_event)
  return (df)
})

ns3_plots <- lapply(SEX_STRATA, function (sx) {
  tmp_df <- dat[[sx]]
  ns3 <- cbind.data.frame(tmp_df$age_event, 
                          ns(tmp_df$scaled_age, 3)) 
  colnames(ns3) <- c("age", "ns1", "ns2", "ns3")
  ns3 <- ns3 %>% pivot_longer(names_to = "ns", values_to = "ns_value", 
                              cols = -age) 
  
  res_plot <- ggplot(ns3, aes(x = age, y = ns_value, 
                              colour = ns)) +
    geom_line() +
    scale_color_brewer(palette = "Set1") +
    labs(x = "Age", y = "Basis value", 
         title = paste0("Natural cubic spline, 3 df for sex: ", sx))
  return (res_plot)
})

# Mixed effects models ----

# Function for best linear/quadratic effect
# Get best model with step-test that performs step-wise elimination
# to determine best random effects and fixed effect structure
get_qmod <- function (df) {
  
  qform <- formula("adj_value ~ 1 + scaled_age + I(scaled_age^2) + (1 + scaled_age +I(scaled_age^2) | eid)")
  q_mod <- tryCatch(lmer(qform, data = df, REML = F),
                    error = function (cond) NA)
  # If the model is too complex, try with simpler random effects
  if (is.na(q_mod)) {
    qform <- formula("adj_value ~ 1 + scaled_age + I(scaled_age^2) + (1 + scaled_age | eid)")
    q_mod <- tryCatch(lmer(qform, data = df, REML = F),
                      error = function (cond) NA)
  }
  # Get best linear + quadratic model
  best_q <- get_model(step(q_mod))
  return (best_q)
}

# Function for spline effect of age
get_smod <- function (df) {
  # Spline effect of age
  sform <- formula("adj_value ~ 1 + splines::ns(scaled_age, 3) + (1 + splines::ns(scaled_age, 3) | eid)")
  s_mod <- tryCatch(lmer(sform, data = df, REML = F),
                    error = function (cond) NA)
  # If the model is too complex, try with simpler random effects
  if (is.na(s_mod)) {
    sform <- formula("adj_value ~ 1 + splines::ns(scaled_age, 3) + (1 | eid)")
    s_mod <- tryCatch(lmer(sform, data = df, REML = F),
                      error = function (cond) NA)
  }
  return (s_mod)
}

all_mods <- lapply(SEX_STRATA, function (sx) {
  # Data for models
  model_df <- dat[[sx]]
  q_mod <- get_qmod(model_df)
  s_mod <- get_smod(model_df)
  
  # Compare and get best model
  to_comp <- c(q_mod, s_mod)
  names(to_comp) <- c("q_mod", "s_mod")
  # Get the best model (BIC)
  mod_BIC <- BIC(q_mod, s_mod)
  best_mod <- to_comp[[rownames(mod_BIC)[which.min(mod_BIC$BIC)]]]
  
  # Save model summary
  sink(log_file, append = T)
  cat("Determining best model structure in sex: ", sx, "\n")
  print(BIC(q_mod, s_mod))
  cat("Best model parameters: ", "\n")
  print(summary(best_mod))
  sink()
  
  # Return best model
  return (best_mod)
  
})
names(all_mods) <- SEX_STRATA

# Save models
saveRDS(all_mods, 
        paste0("best_models_", PHENO, ".rds"))

# Plot residuals ----

qq_plots <- lapply(all_mods, function (bm) {
  resid_dat <- residuals(bm, scaled = T)
  # Down-sample to 10,000 values to plot
  resid_dat <- sample(resid_dat, 10000, replace = F)
  return (qqnorm(resid_dat))
})

prf_plots <- lapply(SEX_STRATA, function (sx) {
  # Residuals vs fitted values
  resid_dat <- data.frame(fitted = fitted(all_mods[[sx]]),
                          resid = residuals(all_mods[[sx]],
                                            scaled = T))
  # Down-sample to 10,000 values to plot
  resid_dat <- resid_dat[sample(1:nrow(resid_dat), 10000, replace = F), ]
  
  prf <- ggplot(resid_dat, aes(x = fitted, y = resid)) +
    geom_point() +
    geom_smooth(method = "lm") +
    labs(x = "fitted value", y = "residual", 
         title = paste0("Residuals plot, sex: ", sx))
  return (prf)
})

# Get terms for clustering and GWAS, and remove outliers ---- 

indiv_dat <- lapply(SEX_STRATA, function (sx) {
  
  res_fixef <- fixef(all_mods[[sx]])
  res_ranef <- ranef(all_mods[[sx]])$eid
  
  # Identify and remove outliers in any term
  outlier_ids <- apply(res_ranef, 2, function (x) {
    outlier_x <- which(x > mean(x) + 5*sd(x) | x < mean(x) - 5*sd(x))
  })
  outlier_ids <- unique(Reduce(union, outlier_ids))
  res_ranef <- res_ranef[-outlier_ids, ]
  # eids to retain
  eids_keep <- rownames(res_ranef)
  
  # Don't keep intercept term for either effect, and only keep fixed effect
  # for portion that needs to be added back into the random term
  RE_NAMES <- colnames(res_ranef)
  RE_NAMES <- RE_NAMES[which(RE_NAMES != "(Intercept)")]
  
  if (length(RE_NAMES) > 0) {
    res_fixef <- res_fixef[names(res_fixef) %in% RE_NAMES]
    # Ensure that the terms align
    res_fixef <- res_fixef[RE_NAMES]
    res_ranef <- as.data.frame(res_ranef[, RE_NAMES])
    colnames(res_ranef) <- RE_NAMES
    
    res_coefs <- apply(res_ranef, 1, function (x) x + res_fixef)
    res_coefs <- as.data.frame(t(res_coefs))
    colnames(res_coefs) <- RE_NAMES
    res_coefs$eid <- eid_order
    
    # Save table of coefficients for clustering and GWAS
    write.table(res_coefs, 
                paste0(PHENO, "/coefficients_", PHENO, "_", sx, ".txt"),
                sep = "\t", col.names = T, quote = F,
                row.names = F)
  } 
  # Return list of IDs that are not outliers
  return (eids_keep)
})

# Get model predictions ----

fitted_dat <- lapply(SEX_STRATA, function (sx) {
  res <- dat[[sx]][, c("eid", "age_event", "scaled_age",
                       "value", "adj_value", "fitted_covs")]
  res$ranef_pred <- predict(all_mods[[sx]])
  res$fixef_pred <- predict(all_mods[[sx]], re.form = NA)
  
  # Add back the covariate effect (which has been adjusted for)
  res <- res %>% 
    mutate(full_ranef = ranef_pred + fitted_covs,
           full_fixef = fixef_pred + fitted_covs)
  
  # Remove outliers (as calculated earlier)
  res <- subset(res, res$eid %in% indiv_dat[[sx]])
  
  write.table(res, 
              paste0("fitted_", PHENO, "_", sx, ".txt"),
              sep = "\t", col.names = T, quote = F)
  return (res)
})
names(fitted_dat) <- SEX_STRATA

# Model plots ----

pdf(paste0("model_plots_", PHENO, ".pdf"), onefile = T)
# Natural cubic spline basis
print(ns3_plots)
# Residuals vs fitted values
print(prf_plots)
# QQ-plots
print(qq_plots)
dev.off()

# Plot a sample of ids to show fixed and random effects ----

pred_pop_level <- lapply(SEX_STRATA, function (sx) {
  rawdat <- dat[[sx]]
  # Remove outliers
  rawdat <- subset(rawdat, rawdat$eid %in% indiv_dat[[sx]])
  
  # Create new data to predict
  newdat <- data.frame(age_event = seq(min(rawdat$age_event), 
                                       max(rawdat$age_event), 
                                       length.out = 100),
                       eid = 0)
  # Get scaled age closest to age from raw data
  get_indexes <- sapply(newdat$age_event, function (x) {
    which.min(abs(x - rawdat$age_event))
  })
  newdat$scaled_age <- rawdat$scaled_age[get_indexes]

  # Population-level prediction with bootstrapped 95% confidence intervals
  pred_pop <- predictInterval(all_mods[[sx]], newdata = newdat,
                              level = 0.95, n.sims = 1000,
                              stat = "median", type = "linear.prediction",
                              include.resid.var = TRUE)
  pred_pop <- cbind(newdat, pred_pop)
  # Add back fixed effect of covariates (median across population)
  COV_EFFECT <- median(rawdat$fitted_covs)
  pred_pop$fit <- pred_pop$fit + COV_EFFECT
  pred_pop$lwr <- pred_pop$lwr + COV_EFFECT
  pred_pop$upr <- pred_pop$upr + COV_EFFECT
  
  # Add sex stratum for combined plot
  pred_pop$sex_strata <- sx
  
  return (pred_pop)
})
names(pred_pop_level) <- SEX_STRATA
pred_pop_level <- bind_rows(pred_pop_level)

# Save models
saveRDS(pred_pop_level, 
        paste0("population_fitted_", PHENO, ".rds"))

# Plot 
pred_plot <- ggplot(pred_pop_level, aes(x = age_event, y = fit,
                                        color = sex_strata,
                                        fill = sex_strata)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) + 
  labs(x = "Age (years)", y = paste(PHENO, "(fitted)"))
