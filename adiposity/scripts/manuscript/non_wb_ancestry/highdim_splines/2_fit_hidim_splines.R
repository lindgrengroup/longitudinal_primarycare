# Author: Samvida S. Venkatesh
# Adapted from: George Nicholson
# Date: 16/05/22

library(splines)
library(lubridate)
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

RANDOM_SEED <- 160522
set.seed(RANDOM_SEED)

plotdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/non_wb_ancestry/highdim_splines/plots/"
resdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/non_wb_ancestry/highdim_splines/"

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

# Load data ----

dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/non_wb_ancestry/highdim_splines/data/dat_to_model.rds")
ANCESTRIES <- names(dat[[1]])

NDF_SPLINE <- 100 # DF of spline
MAX_N_DAYS <- 7500 # Number of days post baseline to be included (~20 years)

# Define spline basis ----

# Create basis
B <- lapply(PHENOTYPES, function (p) {
  per_anc <- lapply(ANCESTRIES, function (anc) {
    per_sex <- lapply(SEX_STRATA, function (sx) {
      diff_day_unique <- 0:min(max(dat[[p]][[anc]][[sx]]$t_diff), MAX_N_DAYS)
      n_days <- length(diff_day_unique)
      res <- splines::bs(diff_day_unique, df = NDF_SPLINE, intercept = TRUE)
      return (res)
    })
    names(per_sex) <- SEX_STRATA
    return (per_sex)
  })
  names(per_anc) <- ANCESTRIES
  return (per_anc)
})
names(B) <- PHENOTYPES

# Modelling functions ----
# All written by GN

#' Function to Calculate AR1 covariance matrix
#' @param n_time Number of time points
#' @param rho AR1 correlation term
#' @param noise_sd Standard deviation of noise process
#' @param intercept_sd Standard deviation of intercept
ar1_covariance <- function(n_time, rho, noise_sd, intercept_sd) {
  R_AR1 <- rho ^ abs(outer(1:n_time, 1:n_time, "-")) # AR1 temporal correlation
  R_diagonal <- diag(n_time)
  R_ones <- matrix(1, n_time, n_time)
  delta_cov <- noise_sd^2 * R_AR1 + intercept_sd^2 * R_ones
  return(delta_cov)
}

#' Calculate subject-specific posterior under smoothing prior
#' @param Z Bayes linear model design matrix
#' @param y Response in linear model
#' @param precision_smooth Precision matrix in MV Gaussian prior on linear model coefficients (zero mean prior currently)
#' @param intercept_sd Standard deviation of intercept
fit_subj_posterior_under_simple_prior <- function(Z, y, precision_smooth) {
  prec_out <- Z %*% t(Z) + precision_smooth
  Sig_out <- solve(prec_out)
  mu_out <- Sig_out %*% Z %*% y
  prec_mu <- prec_out %*% mu_out
  mu_prec_mu <- t(mu_out) %*% prec_mu
  return(list(mu = mu_out, 
              Sig = Sig_out, 
              prec = prec_out, 
              prec_mu = prec_mu, 
              mu_prec_mu = mu_prec_mu,
              log_det_Sig = determinant(Sig_out, log = TRUE)$modulus))
}

# Fit models ----

# Set parameters

AR1_RHO <- 0.99 # Main smoothness parameter
AR1_NOISE_SD <- 2.5 
AR1_INTERCEPT_SD <- 100 # A big number for noninformative intercept in AR1
Sig_smooth <- ar1_covariance(n_time = NDF_SPLINE, 
                             rho = AR1_RHO, 
                             noise_sd = AR1_NOISE_SD, 
                             intercept_sd = AR1_INTERCEPT_SD)
precision_smooth <- solve(Sig_smooth)
str(precision_smooth)

model_dat <- lapply(PHENOTYPES, function (p) {
  per_anc <- lapply(ANCESTRIES, function (anc) {
    per_sex <- lapply(SEX_STRATA, function (sx) {
      sub_dat <- dat[[p]][[anc]][[sx]]
      model_dat <- split(sub_dat, f = sub_dat$eid)
    })
    names(per_sex) <- SEX_STRATA
    return (per_sex)
  })
  names(per_anc) <- ANCESTRIES
  return (per_anc)
})
names(model_dat) <- PHENOTYPES

spline_posteriors <- lapply(PHENOTYPES, function (p) {
  per_anc <- lapply(ANCESTRIES, function (anc) {
    per_sex <- lapply(SEX_STRATA, function (sx) {
      to_return <- lapply(model_dat[[p]][[anc]][[sx]], function (id_df) {
        y <- id_df$value_fulladj_norm
        X <- matrix(0,
                    nrow = min(max(dat[[p]][[anc]][[sx]]$t_diff), MAX_N_DAYS) + 1, 
                    ncol = nrow(id_df))
        for (j in 1:nrow(id_df)) {
          fill_val <- as.numeric(id_df[j, "t_diff"])
          X[fill_val, j] <- 1
        }
        res <- fit_subj_posterior_under_simple_prior(Z = t(B[[p]][[anc]][[sx]]) %*% X,
                                                     y = y,
                                                     precision_smooth = precision_smooth)
        return (res)
      })
      return (to_return)
    })
    names(per_sex) <- SEX_STRATA
    return (per_sex)
  })
  names(per_anc) <- ANCESTRIES
  return (per_anc)
})
names(spline_posteriors) <- PHENOTYPES

# Extract values of interest ----

# To calculate SD_RES, check histogram of residual variances

lapply(PHENOTYPES, function (p) {
  lapply(ANCESTRIES, function (anc) {
    lapply(SEX_STRATA, function (sx) {
      ALL_IDS <- names(model_dat[[p]][[anc]][[sx]])
      resid_vars_check <- lapply(ALL_IDS, function (id) {
        obs_dat <- model_dat[[p]][[anc]][[sx]][[id]]
        pred_dat <- B[[p]][[anc]][[sx]] %*% spline_posteriors[[p]][[anc]][[sx]][[id]]$mu
        
        check <- data.frame(t_diff = obs_dat$t_diff,
                            obs_val = obs_dat$value_fulladj_norm)
        check$pred_val <- pred_dat[check$t_diff]
        var <- 1/nrow(check) * sum((check$obs_val - check$pred_val)^2)
        return (var)
      })
      resid_vars_check <- data.frame(var = unlist(resid_vars_check))
      mean_rvar <- mean(resid_vars_check$var)
      median_rvar <- median(resid_vars_check$var)
      
      resid_var_plot <- ggplot(resid_vars_check, aes(x = var)) +
        geom_histogram() +
        geom_vline(xintercept = mean_rvar, linetype = "dashed") +
        geom_vline(xintercept = median_rvar) +
        labs(x = "Residual variance",
             title = paste0("Mean: ", round(mean_rvar, 3),
                            " Median: ", round(median_rvar, 3)))
      
      pdf(paste0(plotdir, "resid_var_check_", p, "_", anc, "_", sx, ".pdf"),
          onefile = T)
      print(resid_var_plot)
      print(resid_var_plot + scale_x_continuous(limits = c(0, 10)))
      dev.off()
      
      saveRDS(list(B = B[[p]][[anc]][[sx]],
                   spline_posteriors = spline_posteriors[[p]][[anc]][[sx]]),
              paste0(resdir, "fit_objects_", p, "_", anc, "_", sx, ".rds"))
    })
  })
})
