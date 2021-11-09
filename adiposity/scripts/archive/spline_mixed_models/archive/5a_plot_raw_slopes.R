# Author: Samvida S. Venkatesh
# Date: 17/05/21

library(lme4)
library(gamm4)
#library(merTools)
library(tidyverse)
theme_set(theme_bw())
library(RColorBrewer)

set.seed(170521)

# Combine models and covars files across phenotypes and strata ----

PHENO <- c("BMI", "WHR")

## Linear models ----

combined_linear_models <- lapply(PHENO, function (p) {
  # Get saved models
  fnames <- list.files("/well/lindgren/UKBIOBANK/samvida/adiposity/linear/",
                       pattern = paste0(p, "_linear_model_*"))
  # Extract stratum name
  STRATA <- sub(paste0(p, "_linear_model_"), "", fnames)
  STRATA <- sub(".rds", "", STRATA)
  names(fnames) <- STRATA
  # List of models
  models <- lapply(STRATA, function (s) {
    return(readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/linear/",
                          fnames[[s]])))
  })
  names(models) <- STRATA
  return (models)
})
names(combined_linear_models) <- PHENO

saveRDS(combined_linear_models, 
        "/well/lindgren/UKBIOBANK/samvida/adiposity/linear/linear_models.rds")

## Spline models ----

combined_spline_models <- lapply(PHENO, function (p) {
  # Get saved models
  fnames <- list.files("/well/lindgren/UKBIOBANK/samvida/adiposity/splines/",
                       pattern = paste0(p, "_spline_model_*"))
  # Extract stratum name
  STRATA <- sub(paste0(p, "_spline_model_"), "", fnames)
  STRATA <- sub(".rds", "", STRATA)
  names(fnames) <- STRATA
  # List of models
  models <- lapply(STRATA, function (s) {
    return(readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/splines/",
                          fnames[[s]])))
  })
  names(models) <- STRATA
  return (models)
})
names(combined_spline_models) <- PHENO

saveRDS(combined_spline_models, 
        "/well/lindgren/UKBIOBANK/samvida/adiposity/splines/spline_models.rds")

## Slopes files ----

slopes_files <- lapply(PHENO, function (p) {
  # Get saved data
  fnames <- list.files("/well/lindgren/UKBIOBANK/samvida/adiposity/",
                       pattern = paste0(p, "_raw_slopes_and_covars_*"))
  # Extract stratum name
  STRATA <- sub(paste0(p, "_raw_slopes_and_covars_"), "", fnames)
  STRATA <- sub(".rds", "", STRATA)
  names(fnames) <- STRATA
  # List of files
  data <- lapply(STRATA, function (s) {
    return(readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/",
                          fnames[[s]])))
  })
  names(data) <- STRATA
  return (data)
})
names(slopes_files) <- PHENO

saveRDS(slopes_files, 
        "/well/lindgren/UKBIOBANK/samvida/adiposity/raw_slopes_and_covars.rds")

# Read all necessary data ----

lmods <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/linear/linear_models.rds")
smods <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/splines/spline_models.rds")
slope_covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/raw_slopes_and_covars.rds")

# Plot group (fixed effect) predictions (all models) ----

STRATA <- names(lmods[[1]])
PRED_AGES <- seq(20, 80, by = 1)

fe_plots <- lapply(PHENO, function (p) {
  plot_dat <- lapply(STRATA, function (s) {
    
    linear_dat <- data.frame(age_event = PRED_AGES)
    # Get predicted values from linear model (fixed effect)
    linear_dat$mean <- predict(lmods[[p]][[s]], 
                               newdata = linear_dat,
                               re.form = NA)
    # Get 95% CI of fit (from bootstrapping with 100 replications)
    boot_fit <- bootMer(lmods[[p]][[s]], 
                        FUN = function (x) {
                          predict(x, newdata = linear_dat, re.form = NA)}, 
                        nsim = 100, re.form = NA)
    linear_dat$se <- apply(boot_fit$t, 2, sd)
    linear_dat$fit_type <- "linear"
    
    spline_dat <- data.frame(age_event = PRED_AGES)
    # Get predicted values from spline model (fixed effect)
    spline_fit <- predict(smods[[p]][[s]]$gam, 
                          newdata = spline_dat, se.fit = T)
    spline_dat$mean <- spline_fit$fit
    spline_dat$se <- spline_fit$se.fit
    spline_dat$fit_type <- "spline"
    
    res <- bind_rows(linear_dat, spline_dat)
    res$strata <- s
    
    return (res)
  })
  plot_dat <- bind_rows(plot_dat)
  
  plot_res <- ggplot(plot_dat, aes(x = age_event, y = mean,
                                   colour = fit_type, fill = fit_type)) +
    facet_wrap(~strata, ncol = 2, scales = "free") +
    geom_line() +
    geom_ribbon(aes(ymin = mean - 1.96*se,
                    ymax = mean + 1.96*se), alpha = 0.2) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Age (years)", y = "Predicted mean fixed effect (95% CI)",
         title = p)
  
  return (plot_res)
})

pdf("plots/raw_slopes/fixed_effects_model_predictions.pdf", onefile = T)
print(fe_plots)
dev.off()

# Plot sample individual (random effect) fits (all models) ----

re_plots <- lapply(PHENO, function (p) {
  plot_dat <- lapply(STRATA, function (s) {
    
    raw_dat <- lmods[[p]][[s]]@frame
    IDS <- sample(unique(raw_dat$eid), 10, replace = F)
    raw_dat <- subset(raw_dat, raw_dat$eid %in% IDS)
    
    raw_dat$linear_fit <- predict(lmods[[p]][[s]],
                                  newdata = raw_dat)
    
    spline_dat <- smods[[p]][[s]]$mer@frame
    spline_dat <- subset(spline_dat, spline_dat$eid %in% IDS)
    spline_dat$spline_fit <- predict(smods[[p]][[s]]$mer,
                                  newdata = spline_dat)
    
    plot_res <- ggplot(raw_dat, aes(x = age_event, y = value,
                                     group = eid)) +
      geom_point(colour = "#C0C0C0") +
      geom_line(colour = "#C0C0C0") +
      geom_line(aes(y = linear_fit), colour = "#E41A1C") +
      geom_line(data = spline_dat,
                aes(x = age_event, y = spline_fit, group = eid), 
                colour = "#377EB8") +
      labs(x = "Age (years)", y = "Fitted values",
           title = paste(p, s))
    
    return (plot_res)
  })
  return (plot_dat)
})

pdf("plots/raw_slopes/individual_model_predictions.pdf", onefile = T)
print(re_plots)
dev.off()

# Plot correlation between linear and spline BLUPs ----

corr_plots <- lapply(PHENO, function (p) {
  plot_dat <- lapply(STRATA, function (s) {
    res <- slope_covars[[p]][[s]]
    res$strata <- s
    return (res)
  })
  plot_dat <- bind_rows(plot_dat)
  
  plot_res <- ggplot(plot_dat, aes(x = linear_slope, y = spline_slope,
                                   colour = ancestry)) +
    facet_wrap(~strata, ncol = 2, scales = "free") +
    geom_point() +
    scale_colour_brewer(palette = "Dark2", guide = F) + 
    labs(x = "Linear model BLUP", y = "Spline model BLUP",
         title = p)
  
  return (plot_res)
})

pdf("plots/raw_slopes/linear_slope_BLUP_correlation.pdf", onefile = T)
print(corr_plots)
dev.off()

