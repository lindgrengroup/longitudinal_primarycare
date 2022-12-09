# Author: Samvida S. Venkatesh
# Date: 01/11/21

library(argparse)
library(splines)
library(zoo)
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

set.seed(011121)

plots_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/manuscript/figures/"

CLUSTS <- paste0("k", 1:4)

custom_four_diverge <- c("#D35C79", "#D9AB90", "#9FBCA4", "#009593")
names(custom_four_diverge) <- CLUSTS

# Read files ----

parser <- ArgumentParser()
parser$add_argument("--pheno", required=TRUE,
                    help = "Phenotype")
parser$add_argument("--ss", required=TRUE,
                    help = "Sex-strata")
args <- parser$parse_args()

PHENO <- args$pheno
SEX_STRATA <- args$ss

# Define various sets of covariates
TIME_VAR_COVARS <- c("data_provider")
TIME_INVAR_COVARS <- c("baseline_age", "age_sq",
                       "sex", "year_of_birth")

# Centroids from chosen clustering method
clust_centroids <- readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/clustering/",
                                  PHENO, "_", SEX_STRATA, 
                                  "/parameter_selection/K4_L2_M2.rds"))
clust_centroids <- clust_centroids$cluster_centroids

# Soft clustering probabilities
clust_res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/clustering/", 
                                      PHENO, "_", SEX_STRATA,
                               "/soft_clustering_probs_", PHENO, "_", SEX_STRATA, 
                               ".txt"),
                        sep = "\t", header = T, stringsAsFactors = F)

# High-dimensional spline modelling results
hidim_model <- readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/results/with_rvar_fit_objects_", 
                              PHENO, "_", SEX_STRATA, ".rds"))
B <- hidim_model$B
spline_posteriors <- hidim_model$spline_posteriors
model_resid_var <- hidim_model$resid_var

# Covariates
covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[[PHENO]]

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)

# To recapture raw data from residuals 
resid_models <- 
  readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/data/models_for_refitting.rds")[[PHENO]][[SEX_STRATA]]

raw_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/indiv_qcd_data.rds")[[PHENO]]

# Wrangle data -----

covars$eid <- as.character(covars$eid)
general_covars$eid <- as.character(general_covars$eid)

covars <- merge(covars, general_covars, by = "eid")

# Add in covariates to raw GP data
raw_dat$eid <- as.character(raw_dat$eid)
add_covs <- covars %>% 
  select(any_of(c("eid", TIME_VAR_COVARS, TIME_INVAR_COVARS)))
raw_dat <- merge(raw_dat, add_covs, by = "eid")

# Calculate "t" in years
raw_dat <- raw_dat %>% group_by(eid) %>%
  arrange(age_event, .by_group = T) %>%
  mutate(t_diff_yrs = age_event - baseline_age)

MAX_T_DAYS <- 7500
raw_dat <- raw_dat %>% filter(t_diff_yrs <= MAX_T_DAYS/365)

# Plot cluster centroids ----

getPredValuesClusterCentroid <- function (coef_mat) {
  pred_vals <- t(apply(coef_mat, 1, 
                       function (x) B %*% x))
  # Wrangle into ggplot format
  for_plot <- as.data.frame(pred_vals)
  colnames(for_plot) <- paste0("d", 1:ncol(for_plot))
  for_plot$clust <- factor(as.character(1:nrow(for_plot)))
  
  for_plot <- for_plot %>% pivot_longer(cols = -clust,
                                        names_to = "t_diff", 
                                        names_prefix = "d", 
                                        values_to = "pred_value") %>%
    mutate(t_diff = as.numeric(t_diff) / 365)
  return (for_plot)
}

plot_dat <- getPredValuesClusterCentroid(clust_centroids[, -1]) %>%
  mutate(clust = factor(paste0("k", clust),
                        levels = CLUSTS))

clust_centres_plot <- ggplot(plot_dat, 
                             aes(x = t_diff, y = pred_value, 
                                 col = clust, fill = clust)) +
  geom_line() +
  scale_color_manual(values = custom_four_diverge, guide = "none") +
  scale_fill_manual(values = custom_four_diverge, guide = "none") +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE),
                     breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(guide = guide_axis(check.overlap = TRUE),
                     breaks = scales::pretty_breaks(n = 5)) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 6))

tiff(filename = paste0(plots_dir, 
                       PHENO, "_", SEX_STRATA, "_cluster_centroids.tiff"),
     height = 3.75, width = 3.75, units = "cm",
     res = 300)
print(clust_centres_plot)
dev.off()

# Plot sample trajectories ----

ids_per_clust <- lapply(CLUSTS, function (k) {
  thresh_p <- quantile(clust_res[, k], 0.99)
  ids_keep <- clust_res %>% 
    filter(!!as.symbol(k) > thresh_p)
  return (ids_keep$eid)
})
names(ids_per_clust) <- CLUSTS

## Get random IDs within each cluster
getRandIDs <- function (n_each = 5) {
  # Return df of id and cluster number
  sampled_ids <- lapply(CLUSTS, function (k) {
    ret_ids <- sample(ids_per_clust[[k]], n_each, replace = F)
    ret_df <- data.frame(eid = ret_ids,
                         clust = k)
    return (ret_df)
  })
  sampled_ids <- bind_rows(sampled_ids)
  return (sampled_ids)
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

## Create predicted data for set of ids
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
  pred_df$unstd_loci_resid <- 
    ((pred_df$fit_resid - 1.96*pred_df$fit_sd_resid) * sqrt(resid_models$var_full)) + resid_models$mu_full
  pred_df$unstd_upci_resid <- 
    ((pred_df$fit_resid + 1.96*pred_df$fit_sd_resid) * sqrt(resid_models$var_full)) + resid_models$mu_full
  
  
  pred_df <- left_join(pred_df, new_data, by = "eid") %>%
    mutate(t_diff_yrs = (t_diff_days - 1) / 365,
           fit = unstd_fit_resid + add_back_val,
           fit_loci = unstd_loci_resid + add_back_val,
           fit_upci = unstd_upci_resid + add_back_val) %>%
    select(all_of(c("eid", "t_diff_yrs", "fit", "fit_loci", "fit_upci")))
  
  return (pred_df)
}


## Create plots
plotPreds <- function (id_df) {
  sample_raw <- raw_dat %>% filter(eid %in% as.character(id_df$eid))
  sample_pred <- residFitToRaw(as.character(id_df$eid))
  
  sample_pred$clust <- id_df$clust[match(sample_pred$eid,
                                         id_df$eid)]
  sample_pred$clust <- factor(as.character(sample_pred$clust))
  
  # Order the ids by cluster number for plot
  id_levels <- id_df$eid[order(id_df$clust)]
  sample_pred$eid_f <- factor(as.character(sample_pred$eid), 
                              levels = id_levels)
  sample_raw$eid_f <- factor(as.character(sample_raw$eid), 
                             levels = id_levels)
  
  res <- ggplot(sample_pred, aes(x = t_diff_yrs)) +
    facet_wrap(~eid_f, ncol = 4, scales = "free") +
    geom_point(data = sample_raw,
               aes(x = t_diff_yrs, y = value), size = 0.5) +
    geom_line(aes(y = fit,
                  colour = clust)) +
    geom_ribbon(aes(ymin = fit_loci,
                    ymax = fit_upci,
                    fill = clust, colour = clust), alpha = 0.1) +
    scale_color_manual(values = custom_four_diverge, guide = "none") +
    scale_fill_manual(values = custom_four_diverge, guide = "none") +
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

# Samples of trajectories
tiff(filename = paste0(plots_dir, 
                       PHENO, "_", SEX_STRATA, "_clustering_sample_fits.tiff"),
     height = 10, width = 10, units = "cm",
     res = 300)
print(plotPreds(getRandIDs(4)))
dev.off()
