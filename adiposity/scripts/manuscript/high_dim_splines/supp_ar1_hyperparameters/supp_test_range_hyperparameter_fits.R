# Author: Samvida S. Venkatesh
# Adapted from: George Nicholson
# Date: 16/05/22

library(argparse)
library(splines)
library(lubridate)
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

RANDOM_SEED <- 160522
set.seed(RANDOM_SEED)

parser <- ArgumentParser()
parser$add_argument("--phenotype", required = TRUE,
                    help = "Phenotype to model")
parser$add_argument("--sex_strata", required = TRUE,
                    help = "Sex strata")
args <- parser$parse_args()

PHENO <- args$phenotype
SEX_STRATA <- args$sex_strata
NSAMPLES <- 5000

plotdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/ar1_parameter_selection/"
resdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/ar1_parameter_selection/"

dir.create(resdir)
dir.create(plotdir)

# Load data ----

dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/data/dat_to_model_standardised.rds")[[PHENO]][[SEX_STRATA]]

NDF_SPLINE <- 100 # DF of spline
MAX_N_DAYS <- 7500 # Number of days post baseline to be included (~20 years)

# Define spline basis ----

# Create basis
diff_day_unique <- 0:min(max(dat$t_diff), MAX_N_DAYS)
n_days <- length(diff_day_unique)
B <- splines::bs(diff_day_unique, df = NDF_SPLINE, intercept = TRUE)

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

AR1_RHO <- c(0.9, 0.99, 0.999) # Main smoothness parameter
AR1_NOISE_SD <- c(0.5, 2.5, 10) 
AR1_INTERCEPT_SD <- c(10, 100, 500) # A big number for noninformative intercept in AR1

# Create the precision smooths for each combination of parameters
prec_smooths <- lapply(c(1:3), function (set_i) {
  Sig_smooth <- ar1_covariance(n_time = NDF_SPLINE, 
                               rho = AR1_RHO[set_i], 
                               noise_sd = AR1_NOISE_SD[set_i], 
                               intercept_sd = AR1_INTERCEPT_SD[set_i])
  precision_smooth <- solve(Sig_smooth)
  return (precision_smooth)
})

# Select NSAMPLES individuals at random to test sensitivity to hyperparameters
all_ids <- unique(dat$eid)
SAMPLED_IDS <- sample(all_ids, NSAMPLES, replace = F)

dat <- dat %>% filter(eid %in% SAMPLED_IDS)

model_dat <- split(dat, f = dat$eid)

spline_posteriors <- lapply(c(1:3), function (set_i) {
  fitted_dat <- lapply(model_dat, function (id_df) {
    y <- id_df$value_fulladj_norm
    X <- matrix(0, nrow = length(diff_day_unique), ncol = nrow(id_df))
    for (j in 1:nrow(id_df)) {
      fill_val <- as.numeric(id_df[j, "t_diff"])
      X[fill_val, j] <- 1
    }
    res <- fit_subj_posterior_under_simple_prior(Z = t(B) %*% X,
                                                 y = y,
                                                 precision_smooth = prec_smooths[[set_i]])
    return (res)
  })
  return (fitted_dat)
})

# Extract values of interest ----

# To calculate SD_RES, check histogram of residual variances
lapply(c(1:3), function (set_i) {
  resid_vars_check <- lapply(SAMPLED_IDS, function (id) {
    obs_dat <- model_dat[[id]]
    pred_dat <- B %*% spline_posteriors[[set_i]][[id]]$mu
    
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
  
  png(paste0(plotdir, "resid_var_check_parameter_set_", set_i, "_",
             PHENO, "_", SEX_STRATA, ".png"))
  print(resid_var_plot)
  dev.off()
  
  saveRDS(list(B = B,
               spline_posteriors = spline_posteriors[[set_i]],
               resid_var = median_rvar),
          paste0(resdir, "fit_objects_", PHENO, "_", SEX_STRATA, 
                 "_parameter_set_", set_i, ".rds"))
})

