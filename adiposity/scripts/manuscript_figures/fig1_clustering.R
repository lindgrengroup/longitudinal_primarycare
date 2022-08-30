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

# Centroids from chosen clustering method
clust_centroids <- readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/",
                                  PHENO, "_", SEX_STRATA, 
                                  "/parameter_selection/K4_L2_M2.rds"))
clust_centroids <- clust_centroids$cluster_centroids

# Soft clustering probabilities
clust_res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/", 
                                      PHENO, "_", SEX_STRATA,
                               "/soft_clustering_probs_", PHENO, "_", SEX_STRATA, 
                               ".txt"),
                        sep = "\t", header = T, stringsAsFactors = F)

# High-dimensional spline modelling results
hidim_model <- readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/results/fit_objects_", 
                              PHENO, "_", SEX_STRATA, ".rds"))
B <- hidim_model$B
spline_posteriors <- hidim_model$spline_posteriors
model_resid_var <- hidim_model$resid_var

# Covariates
covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[[PHENO]]

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)

# Original data
orig_popn_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/data/dat_to_model.rds")[[PHENO]][[SEX_STRATA]]

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

# Plot population data for individuals w >75% certainty of cluster ----

getSummDat <- function (cluster = "k1", prob_assign = 0.75,
                        valtype = "value", timetype = "age_t1") {
  keep_ids <- as.character(clust_res$eid[which(clust_res[, cluster] >= prob_assign)])
  
  if (timetype == "age_t1") {
    round_yrs = 0.25
    roll_yrs = 2
  } else if (timetype == "t_diff") {
    round_yrs = 100
    roll_yrs = 1000
  }
  
  dat_summ <- orig_popn_dat %>%
    filter(eid %in% keep_ids) %>%
    mutate(time_bin = plyr::round_any(!!as.symbol(timetype), round_yrs, f = round)) %>%
    group_by(time_bin) %>%
    summarise(plot_value = mean_se(!!as.symbol(valtype), 1.96)) %>%
    unnest(plot_value) 
  
  # Get rolling average across 10 years
  dat_summ <- dat_summ %>% 
    mutate(interval_width = seq_along(time_bin) - 
             findInterval(time_bin - roll_yrs, time_bin),
           mean_value_rolled = rollapply(y, interval_width, mean, 
                                         fill = NA),
           lci_value_rolled = rollapply(ymin, interval_width, mean, 
                                        fill = NA),
           uci_value_rolled = rollapply(ymax, interval_width, mean, 
                                        fill = NA))
  to_return <- dat_summ %>%
    select(all_of(c("time_bin", "mean_value_rolled", 
                    "lci_value_rolled", "uci_value_rolled")))
  return (to_return)
}

plotTrajectories <- function (plot_dat) {
  popn_traj_plot <- ggplot(plot_dat,
                           aes(x = time_bin, y = mean_value_rolled, 
                               color = clust, fill = clust)) +
    geom_line() +
    geom_ribbon(aes(ymin = lci_value_rolled, 
                    ymax = uci_value_rolled), linetype = 0,
                alpha = 0.2) +
    scale_color_manual(values = custom_four_diverge, guide = "none") +
    scale_fill_manual(values = custom_four_diverge, guide = "none") +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE),
                       breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(guide = guide_axis(check.overlap = TRUE),
                       breaks = scales::pretty_breaks(n = 5)) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_text(size = 6))
  return (popn_traj_plot)
}

# Apply to adj- and baselined- data

to_plot_adj_baselined <- lapply(CLUSTS, function (clustk) {
  summ_clust <- getSummDat(cluster = clustk, 
                           prob_assign = 0.75,
                           valtype = "value_fulladj", 
                           timetype = "t_diff") %>%
    mutate(clust = clustk)
  
  return (summ_clust)
})
to_plot_adj_baselined <- bind_rows(to_plot_adj_baselined) %>%
  mutate(clust = factor(clust, 
                        levels = CLUSTS))

# To convert t_diff to yrs
to_plot_adj_baselined <- to_plot_adj_baselined %>%
  mutate(time_bin = as.numeric(time_bin) / 365)

tiff(filename = paste0(plots_dir, 
                       PHENO, "_", SEX_STRATA, "_popn_clusters.tiff"),
     height = 3.75, width = 3.75, units = "cm",
     res = 300)
print(plotTrajectories(to_plot_adj_baselined))
dev.off()

# Plot sample trajectories ----

# Assign individuals to a single cluster 
conf_clusts <- clust_res %>%
  filter(if_any(starts_with("k"), any_vars(. >= 0.75)))

clust_assts <- conf_clusts %>%
  select(-eid) %>%
  mutate(clust = names(.)[max.col(.)]) 
clust_assts$eid <- conf_clusts$eid

# Count number of individuals in each cluster to make sampling quicker later
clust_assts <- clust_assts %>% 
  group_by(clust) %>% 
  mutate(n_cluster = n())

## Get random IDs within each cluster
getRandIDs <- function (n_each = 5) {
  # Return df of id and cluster number
  sampled_ids <- clust_assts %>% 
    group_by(clust) %>%
    sample_n(ifelse(n_cluster < n_each, n_cluster, n_each)) %>%
    mutate(eid = as.character(eid)) %>%
    select(eid, clust)
  return (sampled_ids)
}

## Create predicted data for set of ids
createPredDat <- function (id_list) {
  res_list <- lapply(id_list, function (id) {
    pred_value <- B %*% spline_posteriors[[id]]$mu
    sd_pred <- sqrt(diag(B %*% spline_posteriors[[id]]$Sig %*% t(B)) * model_resid_var)
    
    res <- data.frame(eid = id,
                      t_diff = 1:length(pred_value),
                      fit_mean = pred_value,
                      fit_sd = sd_pred)
  })
  res_list <- bind_rows(res_list)
  return (res_list)
}

## Create plots
plotPreds <- function (id_df) {
  sample_raw <- orig_popn_dat %>% filter(eid %in% id_df$eid)
  sample_pred <- createPredDat(id_df$eid)
  
  # Make t_diff into yrs
  sample_raw <- sample_raw %>%
    mutate(t_diff_yrs = (t_diff - 1)/365)
  sample_pred <- sample_pred %>%
    mutate(t_diff_yrs = (t_diff - 1)/365) 
  
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
    facet_wrap(~eid_f, ncol = 4, scales = "free_y") +
    geom_point(data = sample_raw,
               aes(x = t_diff_yrs, y = value_fulladj), size = 0.5) +
    geom_line(aes(y = fit_mean,
                  colour = clust)) +
    geom_ribbon(aes(ymin = fit_mean - 1.96*fit_sd,
                    ymax = fit_mean + 1.96*fit_sd,
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
