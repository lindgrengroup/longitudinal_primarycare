# Author: Samvida S. Venkatesh
# Date: 01/11/21

library(argparse)
library(lme4)
library(splines)
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

set.seed(011121)

infile_path <- "" # REDACTED
gen_resources_path <- "" # REDACTED
plots_dir <- "" # REDACTED

# For two colours (rose and navy):
custom_two_diverge <- c("#D35C79", "#005580")

# Read files ----

parser <- ArgumentParser()
parser$add_argument("--pheno", required=TRUE,
                    help = "Phenotype")
parser$add_argument("--ss", required=TRUE,
                    help = "Sex-strata")
args <- parser$parse_args()

PHENO <- args$pheno
SEX_STRATA <- args$ss

lmm_model <- readRDS(paste0(infile_path, "/lmm_models/",
                            PHENO, "_full_model.rds"))[[SEX_STRATA]]

# High-dimensional spline modelling results
hidim_model <- readRDS(paste0(infile_path, "/highdim_splines/standardised_outcomes/results/with_rvar_fit_objects_", 
                              PHENO, "_", SEX_STRATA, ".rds"))
B <- hidim_model$B
spline_posteriors <- hidim_model$spline_posteriors
model_resid_var <- hidim_model$resid_var

# To recapture raw data from residuals 
resid_models <- 
  readRDS(psate0(infile_path, "/highdim_splines/data/models_for_refitting.rds"))[[PHENO]][[SEX_STRATA]]

raw_dat <- readRDS(paste0(infile_path, "/data/indiv_qcd_data.rds"))[[PHENO]]

# Covariates
covars <- readRDS(paste0(infile_path, "/data/covariates.rds"))[[PHENO]]

general_covars <- read.table(paste0(gen_resources_path, "/220504_QCd_demographic_covariates.txt"),
                             sep = "\t", header = T, stringsAsFactors = F)

# Define various sets of covariates
TIME_VAR_COVARS <- c("data_provider")
TIME_INVAR_COVARS <- c("baseline_age", "age_sq",
                       "sex", "year_of_birth")

# Wrangle data ----

covars$eid <- as.character(covars$eid)
general_covars$eid <- as.character(general_covars$eid)

covars <- merge(covars, general_covars, by = "eid")

# Add in covariates to raw GP data
raw_dat$eid <- as.character(raw_dat$eid)
add_covs <- covars %>% 
  select(any_of(c("eid", TIME_VAR_COVARS, TIME_INVAR_COVARS)))
raw_dat <- merge(raw_dat, add_covs, by = "eid")

# Calculate "t" in years
raw_dat <- raw_dat %>% mutate(t = age_event - baseline_age)

# Sample individuals on whom to run all fitting procedures ----

all_ids_hidim <- as.character(unique(raw_dat$eid))
all_ids_lmm <- as.character(unique(lmm_model@frame$eid))
ALL_IDS <- intersect(all_ids_hidim, all_ids_lmm)

RETAIN_IDS <- sample(ALL_IDS, 16, replace = F)
raw_dat <- raw_dat %>% filter(eid %in% RETAIN_IDS)

# How far out to predict (in days)
MAX_T_DAYS <- 7500
raw_dat <- raw_dat %>% filter(t <= MAX_T_DAYS/365)
spline_posteriors <- spline_posteriors[RETAIN_IDS]

# Fitting functions ----

fitLMM <- function (id_list) {
  
  # Create new data to predict from 
  new_data <- raw_dat %>%
    filter(eid %in% id_list) %>%
    # Get covariates 
    select(any_of(c("eid", TIME_VAR_COVARS, TIME_INVAR_COVARS))) %>% 
    distinct(across(all_of(c("eid", TIME_VAR_COVARS))), .keep_all = T) 
  
  # Timepoints to extend to
  # Predict for every 1 yr interval up to the max num-yrs
  timepoints <- seq(0, MAX_T_DAYS/365, by = 1)
  nts <- nrow(new_data)
  new_data <- new_data %>% slice(rep(1:n(), each = length(timepoints)))
  
  new_data <- new_data %>% 
    mutate(t = rep(timepoints, times = nts),
           age_event = t + baseline_age,
           age_event_sq = age_event^2)
  
  # Predict new values
  fitted_results <- as.data.frame(predict(lmm_model,
                                          newdata = new_data))
  colnames(fitted_results) <- "fit"
  pred_df <- bind_cols(new_data, fitted_results)
  
  # At each time-point, average across time-varying covars
  pred_df <- pred_df %>% 
    group_by(eid, t) %>%
    summarise(fit = mean(fit),
              fit_sd = NA) %>%
    rename(t_diff_yrs = t) %>% ungroup() %>%
    select(all_of(c("eid", "t_diff_yrs", "fit", "fit_sd")))
  
  return (pred_df)
}

fitHiDimSpline <- function (id_list) {
  res_list <- lapply(id_list, function (id) {
    pred_value <- B %*% spline_posteriors[[id]]$mu
    sd_pred <- sqrt(diag(B %*% spline_posteriors[[id]]$Sig %*% t(B)) * model_resid_var)
    
    res <- data.frame(eid = id,
                      t_diff_days = 1:length(pred_value),
                      fit_resid = pred_value,
                      fit_sd_resid = sd_pred)
    res <- res %>% filter(t_diff_days <= MAX_T_DAYS)
    return (res)
  })
  res_list <- bind_rows(res_list)
  return (res_list)
}

# Convert high-dimensional spline residual back to original value by 
# un-scaling and adding predicted

residFitToRaw <- function (id_list) {
  # Create new data to predict from 
  new_data <- raw_dat %>%
    filter(eid %in% id_list) %>%
    # Get covariates 
    select(any_of(c("eid", TIME_VAR_COVARS, TIME_INVAR_COVARS))) %>% 
    distinct(across(all_of(c("eid", TIME_VAR_COVARS))), .keep_all = T) 
  # Get predicted value (to add back to fitted residual, averaging over time-varying covars)
  new_data$add_back_val <- predict(resid_models$fullmod, 
                                   newdata = new_data)
  new_data <- new_data %>% 
    group_by(eid) %>%
    summarise(add_back_val = mean(add_back_val, na.rm = T))
  
  # Predict new values
  pred_df <- fitHiDimSpline(id_list)
  # Un-standardise fitted values
  pred_df$unstd_fit_resid <- 
    (pred_df$fit_resid * sqrt(resid_models$var_full)) + resid_models$mu_full
  
  pred_df <- left_join(pred_df, new_data, by = "eid") %>%
    mutate(t_diff_yrs = (t_diff_days - 1) / 365,
           fit = unstd_fit_resid + add_back_val,
           fit_sd = unstd_fit_resid + add_back_val) %>%
    select(all_of(c("eid", "t_diff_yrs", "fit", "fit_sd")))
  
  return (pred_df)
}

# Plotting functions ----

# Given dataframe with fits and fit types, make plots for predictions
plotPredictions <- function (raw_df, lmm_fits, hidim_fits) {
  
  raw_df <- raw_df %>% rename(t_diff_yrs = t)
  lmm_fits$model_type <- "lmm"
  hidim_fits$model_type <- "hidim_spline"
  fit_df <- bind_rows(lmm_fits, hidim_fits)
  
  # Get max time to plot for each individual
  max_t_dat <- raw_df %>% group_by(eid) %>% 
    summarise(max_t = max(t_diff_yrs))
  fit_df <- fit_df %>% 
    mutate(max_t_plot = max_t_dat$max_t[match(eid, max_t_dat$eid)]) %>%
    filter(t_diff_yrs <= max_t_plot)
  
  res <- ggplot(raw_df, aes(x = t_diff_yrs)) +
    facet_wrap(~new_id, ncol = 4, scales = "free") +
    # Plot raw data underneath
    geom_point(data = raw_df, 
               aes(y = value),
               colour = "#404040", size = 0.7) +
    # Ribbon for 95% CI
    geom_line(data = fit_df, aes(y = fit, colour = model_type)) +
    # geom_ribbon(data = fit_df, aes(ymin = fit - fit_sd,
    #                 ymax = fit + fit_sd,
    #                 colour = model_type, fill = model_type), 
    #             alpha = 0.1, linetype = 0) +
    scale_color_manual(values = custom_two_diverge, guide = F) +
    scale_fill_manual(values = custom_two_diverge, guide = F) +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE),
                       breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(guide = guide_axis(check.overlap = TRUE),
                       breaks = scales::pretty_breaks(n = 3)) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_text(size = 6),
          strip.background = element_blank(),
          strip.text.x = element_blank())
  
  return (res)
}

# Apply plotting functions ----

lmm_fits <- fitLMM(RETAIN_IDS)
hidim_fits <- residFitToRaw(RETAIN_IDS)

# Rename all ids for anonymity
old_to_new <- data.frame(old_id = as.character(RETAIN_IDS),
                         new_id = paste0("id_", 1:length(RETAIN_IDS)))

plot_raw_dat <- raw_dat %>% 
  mutate(new_id = old_to_new$new_id[match(eid, old_to_new$old_id)],
         new_id = factor(new_id, levels = paste0("id_", 1:length(RETAIN_IDS))))
plot_lmm_fits <- lmm_fits %>% 
  mutate(new_id = old_to_new$new_id[match(eid, old_to_new$old_id)],
         new_id = factor(new_id, levels = paste0("id_", 1:length(RETAIN_IDS))))
plot_hidim_fits <- hidim_fits %>% 
  mutate(new_id = old_to_new$new_id[match(eid, old_to_new$old_id)],
         new_id = factor(new_id, levels = paste0("id_", 1:length(RETAIN_IDS))))

tiff(filename = paste0(plots_dir, 
                       PHENO, "_", SEX_STRATA, "_sample_fits.tiff"),
     height = 10, width = 10, units = "cm",
     res = 300)
print(plotPredictions(raw_df = plot_raw_dat,
                      lmm_fits = plot_lmm_fits,
                      hidim_fits = plot_hidim_fits))
dev.off()
