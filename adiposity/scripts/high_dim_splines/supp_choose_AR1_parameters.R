# Author: Samvida S. Venkatesh
# Adapted from: George Nicholson
# Date: 14/11/22

library(splines)
library(lubridate)
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

RANDOM_SEED <- 160522
set.seed(RANDOM_SEED)

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")
plotdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/ar1_parameter_selection/"

# Load data ----

dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/data/dat_to_model_standardised.rds")
covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")

NDF_SPLINE <- 100 # DF of spline
MAX_N_DAYS <- 7500 # Number of days post baseline to be included (~20 years)

# Define spline basis ----

# Create basis
B <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    diff_day_unique <- 0:min(max(dat[[p]][[sx]]$t_diff), MAX_N_DAYS)
    n_days <- length(diff_day_unique)
    basis <- splines::bs(diff_day_unique, df = NDF_SPLINE, intercept = TRUE)
    return (basis)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(B) <- PHENOTYPES

# Select subset of IDs to visualise modelling results ----

NSAMPLE <- 25
NQS <- 5 # number of quantiles to sample within

ids_select <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    # Select NSAMPLE/NQS ids from each quantile of # FU measures
    get_ids <- covars[[p]] %>% filter(eid %in% dat[[p]][[sx]]$eid) %>%
      mutate(FU_n_gp = cut_number(FU_n, NQS)) %>%
      group_by(FU_n_gp) %>%
      slice_sample(n = NSAMPLE/NQS) %>%
      select(eid, FU_n_gp)
    
    # Reassign ids based on group
    get_ids <- get_ids %>%
      arrange(FU_n_gp) %>% ungroup() %>%
      mutate(new_id = factor(paste0("id_", row_number()),
                             levels = paste0("id_", 1:NSAMPLE)))
    return (get_ids)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
}) 
names(ids_select) <- PHENOTYPES

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

# Functions to visualise modelling results ----
# Without standard deviation for now, because this will need to be calculated
# from residual variance across all individuals

getPredDat <- function (p, sx, spline_posts) {
  ids_fit <- ids_select[[p]][[sx]]
  B_basis <- B[[p]][[sx]]
  
  pred_sample <- lapply(ids_fit$eid, function (id) {
    pred_value <- B_basis %*% spline_posts[[id]]$mu
    res <- data.frame(eid = id,
                      t_diff = 1:length(pred_value),
                      fit_mean = pred_value)
    return (res)
  })
  pred_sample <- bind_rows(pred_sample)
  pred_sample$plot_id <- ids_fit$new_id[match(pred_sample$eid,
                                              ids_fit$eid)]
  return (pred_sample)
}

plotPreds <- function (p, sx, pdat, varying_par) {
  ids_fit <- ids_select[[p]][[sx]]
  raw_sample <- dat[[p]][[sx]] %>% filter(eid %in% ids_fit$eid)
  raw_sample$plot_id <- ids_fit$new_id[match(raw_sample$eid,
                                             ids_fit$eid)]
  
  if (varying_par == "ar1_rho") {
    colpal_use <- brewer.pal(n = 9, "Blues")[3:9]
  } else if (varying_par == "ar1_noise_sd") {
    colpal_use <- brewer.pal(n = 9, "Oranges")[3:9] 
  } else if (varying_par == "ar1_isd") { 
    colpal_use <- brewer.pal(n = 9, "Greens")[3:9]
  } else {
    colpal_use <- brewer.pal(n = 9, "Greys")[3:9]
  }
  
  # Plot
  if (varying_par != "none") {
    pdat <- pdat %>%
      mutate(vpar = factor(!!as.symbol(varying_par)))
    fit_plot <- ggplot(pdat, aes(x = t_diff)) +
      facet_wrap(~plot_id, nrow = NQS, scales = "free_y") +
      geom_point(data = raw_sample,
                 aes(x = t_diff, y = value_fulladj_norm)) +
      geom_line(aes(y = fit_mean, colour = vpar)) +
      scale_color_manual(values = colpal_use) + 
      scale_x_continuous(guide = guide_axis(check.overlap = TRUE),
                         breaks = scales::pretty_breaks(n = 3)) +
      scale_y_continuous(guide = guide_axis(check.overlap = TRUE),
                         breaks = scales::pretty_breaks(n = 3)) +
      theme(axis.text = element_text(size = 6)) +
      labs(x = "Days from first measurement", y = "Adj. standardised value",
           title = paste0("Varying ", varying_par))
  } else {
    fit_plot <- ggplot(pdat, aes(x = t_diff)) +
      facet_wrap(~plot_id, nrow = NQS) +
      geom_point(data = raw_sample,
                 aes(x = t_diff, y = value_fulladj_norm)) +
      geom_line(aes(y = fit_mean)) +
      scale_x_continuous(guide = guide_axis(check.overlap = TRUE),
                         breaks = scales::pretty_breaks(n = 3)) +
      scale_y_continuous(limits = c(-1.5, 1.5),
                         breaks = scales::pretty_breaks(n = 3)) +
      theme(axis.text = element_text(size = 6)) +
      labs(x = "Days from first measurement", y = "Adj. standardised value",
           title = paste0("Varying ", varying_par))
  }
  
  return (fit_plot)
}

# Fit models ----

# Parameters to try
AR1_RHO <- c(0.99, 0.9925, 0.995, 0.975, 0.999) # Main smoothness parameter
AR1_NOISE_SD <- c(2, 5, 10, 20) # Overfitting vs underfitting control
AR1_INTERCEPT_SD <- 100 # A big number for noninformative intercept in AR1

# Create the precision smooths for each combination of parameters
prec_smooths <- lapply(AR1_RHO, function (ar1r) {
  per_sigsq <- lapply(AR1_NOISE_SD, function (sigsq) {
    per_isd <- lapply(AR1_INTERCEPT_SD, function (ar1_isd) {
      Sig_smooth <- ar1_covariance(n_time = NDF_SPLINE, 
                                   rho = ar1r, 
                                   noise_sd = sigsq, 
                                   intercept_sd = ar1_isd)
      precision_smooth <- solve(Sig_smooth)
      return (precision_smooth)
    })
    names(per_isd) <- paste0("ar1_isd_", AR1_INTERCEPT_SD)
    return (per_isd)
  })
  names(per_sigsq) <- paste0("sigsq_", AR1_NOISE_SD)
  return (per_sigsq)
})
names(prec_smooths) <- paste0("ar1_rho_", AR1_RHO)

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    # Get sample of ids to test modelling on
    sub_dat <- dat[[p]][[sx]] %>%
      filter(eid %in% ids_select[[p]][[sx]]$eid)
    model_dat <- split(sub_dat, f = sub_dat$eid)
    all_ids <- names(model_dat)
    B_basis <- B[[p]][[sx]] 
    
    full_res <- lapply(AR1_RHO, function (ar1r) {
      per_sigsq <- lapply(AR1_NOISE_SD, function (sigsq) {
        per_isd <- lapply(AR1_INTERCEPT_SD, function (ar1_isd) {
          
          prec_smooth_use <- prec_smooths[[paste0("ar1_rho_", ar1r)]][[paste0("sigsq_", sigsq)]][[paste0("ar1_isd_", ar1_isd)]]
          
          spline_posteriors <- lapply(model_dat, function (id_df) {
            y <- id_df$value_fulladj_norm
            X <- matrix(0, nrow = nrow(B_basis), ncol = nrow(id_df))
            for (j in 1:nrow(id_df)) {
              fill_val <- as.numeric(id_df[j, "t_diff"])
              X[fill_val, j] <- 1
            }
            res <- fit_subj_posterior_under_simple_prior(Z = t(B_basis) %*% X,
                                                         y = y,
                                                         precision_smooth = prec_smooth_use)
            return (res)
          })
          pdat <- getPredDat(p, sx, spline_posteriors)
          pdat$ar1_isd <- ar1_isd
          return (pdat)
        })
        per_isd <- bind_rows(per_isd)
        per_isd$ar1_noise_sd <- sigsq
        return (per_isd)
      })
      per_sigsq <- bind_rows(per_sigsq)
      per_sigsq$ar1_rho <- ar1r
      return (per_sigsq)
    })
    full_res <- bind_rows(full_res)
    
    # Plot results coloured by various aspects
    
    # Fix NOISE_SD
    pdat <- full_res %>% filter(ar1_noise_sd == 5)
    png(paste0(plotdir, p, "_", sx, "_fits_sigsq_5_ar1_isd_100.png"),
        res = 300, units = "cm", height = 20, width = 25)
    print(plotPreds(p, sx, pdat, varying_par = "ar1_rho"))
    dev.off()
    
    # Fix RHO
    pdat <- full_res %>% filter(ar1_rho == 0.995)
    png(paste0(plotdir, p, "_", sx, "_fits_ar1_rho_0.995_ar1_isd_100.png"),
        res = 300, units = "cm", height = 20, width = 25)
    print(plotPreds(p, sx, pdat, varying_par = "ar1_noise_sd"))
    dev.off()
    
    # FIX ALL
    pdat <- full_res %>% filter(ar1_rho == 0.995 & ar1_noise_sd == 10)
    png(paste0(plotdir, p, "_", sx, "_fits_ar1_rho_0.995_sigsq_10_ar1_isd_100.png"),
        res = 300, units = "cm", height = 20, width = 25)
    print(plotPreds(p, sx, pdat, varying_par = "none"))
    dev.off()
    
  })
})
