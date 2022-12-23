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

mainpath <- "" # REDACTED
gen_resources_path <- "" # REDACTED
hidim_mods_path <- "" # REDACTED

plotdir <- paste0(hidim_mods_path, "/plots/")
resdir <- paste0(hidim_mods_path, "/results/")

dir.create(plotdir)
dir.create(resdir)

# Load data ----

dat <- readRDS(paste0(hidim_mods_path, "/data/dat_to_model_standardised.rds"))[[PHENO]][[SEX_STRATA]]

NDF_SPLINE <- 100 # DF of spline
MAX_N_DAYS <- 7500 # Number of days post baseline to be included (~20 years)

# Sanity check - plot sample of data for 25 individuals ----

# set.seed(RANDOM_SEED)
# plot_indivs <- sample(unique(dat$eid), 25, replace = F)
# plot_dat <- dat %>% filter(eid %in% plot_indivs)
# 
# rawdat_plot <- ggplot(plot_dat, aes(x = t_diff, y = value)) +
#   facet_wrap(~eid, nrow = 5, ncol = 5, scales = "free") +
#   geom_point() +
#   labs(x = "Days from first measurement", y = "Observed value")
# 
# minadj_plot <- ggplot(plot_dat, aes(x = t_diff, y = value_minadj)) +
#   facet_wrap(~eid, nrow = 5, ncol = 5, scales = "free") +
#   geom_point() +
#   labs(x = "Days from first measurement", y = "Age- and sex-adj value")
# 
# fulladj_plot <- ggplot(plot_dat, aes(x = t_diff, y = value_fulladj)) +
#   facet_wrap(~eid, nrow = 5, ncol = 5, scales = "free") +
#   geom_point() +
#   labs(x = "Days from first measurement", y = "Confounder-adj value")
# 
# pdf(paste0(plotdir, "sample_obs_", PHENO, "_", SEX_STRATA, ".pdf"), onefile = T)
# print(rawdat_plot)
# print(minadj_plot)
# print(fulladj_plot)
# dev.off()

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

AR1_RHO <- 0.99 # Main smoothness parameter
AR1_NOISE_SD <- 2.5 
AR1_INTERCEPT_SD <- 100 # A big number for noninformative intercept in AR1
Sig_smooth <- ar1_covariance(n_time = NDF_SPLINE, 
                             rho = AR1_RHO, 
                             noise_sd = AR1_NOISE_SD, 
                             intercept_sd = AR1_INTERCEPT_SD)
precision_smooth <- solve(Sig_smooth)
str(precision_smooth)

model_dat <- split(dat, f = dat$eid)
ALL_IDS <- names(model_dat)

spline_posteriors <- lapply(model_dat, function (id_df) {
  y <- id_df$value_fulladj_norm
  X <- matrix(0, nrow = length(diff_day_unique), ncol = nrow(id_df))
  for (j in 1:nrow(id_df)) {
    fill_val <- as.numeric(id_df[j, "t_diff"])
    X[fill_val, j] <- 1
  }
  res <- fit_subj_posterior_under_simple_prior(Z = t(B) %*% X,
                                               y = y,
                                               precision_smooth = precision_smooth)
  return (res)
})

# Extract values of interest ----

# To calculate SD_RES, check histogram of residual variances
resid_vars_check <- lapply(ALL_IDS, function (id) {
  obs_dat <- model_dat[[id]]
  pred_dat <- B %*% spline_posteriors[[id]]$mu
  
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

png(paste0(plotdir, "resid_var_check_", PHENO, "_", SEX_STRATA, ".png"))
print(resid_var_plot)
dev.off()

saveRDS(list(B = B,
             spline_posteriors = spline_posteriors,
             resid_var = median_rvar),
        paste0(resdir, "fit_objects_", PHENO, "_", SEX_STRATA, ".rds"))

# Plot fitted values for sample data ----

## only run this after getting the residual variance from above

# set.seed(RANDOM_SEED)
# plot_indivs <- sample(unique(dat$eid), 25, replace = F)
# 
# pred_sample <- lapply(plot_indivs, function (id) {
#   pred_value <- B %*% spline_posteriors[[id]]$mu
#   sd_pred <- sqrt(diag(B %*% spline_posteriors[[id]]$Sig %*% t(B)) * VAR_RES)
#   
#   res <- data.frame(eid = id,
#                     t_diff = 1:length(pred_value),
#                     fit_mean = pred_value,
#                     fit_sd = sd_pred)
#   return (res)
# })
# pred_sample <- bind_rows(pred_sample)
# 
# raw_sample <- dat %>% filter(eid %in% plot_indivs)
# 
# fit_plot <- ggplot(pred_sample, aes(x = t_diff)) +
#   facet_wrap(~eid, nrow = 5, ncol = 5, scales = "free_y") +
#   geom_point(data = raw_sample,
#              aes(x = t_diff, y = value_fulladj)) +
#   geom_line(aes(y = fit_mean)) +
#   geom_ribbon(aes(ymin = fit_mean - 1.96*fit_sd,
#                   ymax = fit_mean + 1.96*fit_sd), alpha = 0.1) +
#   labs(x = "Days from first measurement", y = "Confounder-adj value")
# 
# pdf(paste0(plotdir, "sample_pred_", PHENO, "_", SEX_STRATA, ".pdf"), onefile = T)
# print(fit_plot)
# dev.off()

# # Matrix of means (id x basis) 
# mn_mat <- lapply(spline_posteriors, function (spobj) {
#   return (as.data.frame(t(spobj$mu)))
# })
# mn_mat <- bind_rows(mn_mat)
# rownames(mn_mat) <- names(spline_posteriors)
# 
# # Matrix of standard deviations (id x basis)
# sd_mat <- lapply(spline_posteriors, function (spobj) {
#   return (sqrt(diag(spobj$Sig * VAR_RES)))
# })
# sd_mat <- bind_rows(sd_mat)
# sd_mat <- as.data.frame(sd_mat)
# rownames(sd_mat) <- names(spline_posteriors)
# # Original data with fitted values and SDs
# fitted_dat <- lapply(ALL_IDS, function (id) {
#   res <- model_dat[[id]]
#   pred_value <- B %*% spline_posteriors[[id]]$mu
#   sd_pred <- sqrt(diag(B %*% spline_posteriors[[id]]$Sig %*% t(B)) * VAR_RES)
#   
#   res$fitted_mean_value_fulladj <- pred_value[res$t_diff]
#   res$fitted_sd_value_fulladj <- sd_pred[res$t_diff]
#   return (res)
# })
# fitted_dat <- bind_rows(fitted_dat)