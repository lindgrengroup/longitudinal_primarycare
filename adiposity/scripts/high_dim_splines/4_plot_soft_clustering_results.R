# Author: Samvida S. Venkatesh
# Date: 24/05/2022

library(argparse)
library(splines)
library(tidyverse)
library(zoo)
library(pheatmap)
library(ggpubr)
library(GGally)
library(RColorBrewer)
theme_set(theme_bw())

RANDOM_SEED <- 160522
set.seed(RANDOM_SEED)
custom_four_diverge <- c("#D35C79", "#D9AB90", "#9FBCA4", "#009593")
names(custom_four_diverge) <- c("k1", "k2", "k3", "k4")

# Get arguments ----

parser <- ArgumentParser()
parser$add_argument("--phenotype", required = TRUE,
                    help = "Phenotype to model")
parser$add_argument("--sex_strata", required = TRUE,
                    help = "Sex strata")

args <- parser$parse_args()

PHENO <- args$phenotype
SEX_STRATA <- args$sex_strata
K_chosen <- 4

plotdir <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/", 
                  PHENO, "_", SEX_STRATA, "/plots/soft_clustering/")
dir.create(plotdir)
resdir <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/", 
                 PHENO, "_", SEX_STRATA, "/")

# Load data ----

# High-dimensional spline modelling results
model_dat <- readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/results/fit_objects_", 
                            PHENO, "_", SEX_STRATA, ".rds"))
B <- model_dat$B
spline_posteriors <- model_dat$spline_posteriors
model_resid_var <- model_dat$resid_var

# Soft clustering probabilities
clust_res <- read.table(paste0(resdir, "soft_clustering_probs_", PHENO, "_", SEX_STRATA, 
                               ".txt"),
                        sep = "\t", header = T, stringsAsFactors = F)

# Covariates
covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[[PHENO]]

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
covars$eid <- as.character(covars$eid)
general_covars$eid <- as.character(general_covars$eid)

covars <- merge(covars, general_covars, by = "eid")

# Original data
orig_popn_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/data/dat_to_model.rds")[[PHENO]][[SEX_STRATA]]

# Plot soft clustering results ----

CAT_COVARS <- "sex"
QUANT_COVARS <- c("baseline_age", "baseline_trait",
                  "FUyrs", "FU_n",
                  "year_of_birth", "age_at_death")

## Plots for max probability per individual -----

# Histogram
for_hist <- clust_res %>%
  mutate(max_prob = dplyr::select(., starts_with("k")) %>% 
           do.call(pmax, .))

hist_plot <- ggplot(for_hist, aes(x = max_prob)) +
  geom_histogram() +
  expand_limits(x = 0)

png(paste0(plotdir, "max_probabilities_hist.png"))
print(hist_plot)
dev.off()

# vs covariates
for_hist <- for_hist %>%
  dplyr::select(all_of(c("eid", "max_prob"))) %>%
  # Cut max prob into 0.1 unit bins (<0.5, <0.6, etc.)
  mutate(max_prob = as.factor(as.character(paste0("<", 
                                                  plyr::round_any(max_prob, 0.1, f = floor)))),
         eid = as.character(eid))

covar_hists <- left_join(for_hist, covars, by = "eid")
plot_list <- lapply(QUANT_COVARS, function (cov) {
  res <- ggplot(covar_hists %>% filter(!is.na(!!as.symbol(cov))), 
                aes(x = max_prob, y = !!as.symbol(cov))) +
    geom_violin(aes(fill = max_prob), 
                position = position_dodge(1)) +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "max prob of k", y = cov) +
    theme(legend.position = "none")
  return (res)
})
arranged_plots <- ggarrange(plotlist = plot_list, nrow = 3, ncol = 2)
png(paste0(plotdir, "max_probabilities_vs_covars.png"))
print(arranged_plots)
dev.off()

## Plots for probability of one cluster vs another -----

pairs_plot <- ggpairs(clust_res, columns = 2:(K_chosen+1))
png(paste0(plotdir, "clust_probs_pairs_plot.png"))
print(pairs_plot)
dev.off()

## Heatmaps ----

# B-coefficients mean matrix
mn_mat <- lapply(spline_posteriors, function (spobj) {
  return (as.data.frame(t(spobj$mu)))
})
mn_mat <- bind_rows(mn_mat)
rownames(mn_mat) <- names(spline_posteriors)

# Function to RINT coefficient for plotting
RINTx <- function (x) {
  return (qnorm((rank(x) - 0.5) / sum(!is.na(x))))
}

# Function to order ids by probability of belonging to cluster "k"
orderIDs <- function (cluster = "k1") {
  return (order(clust_res[, cluster], decreasing = T))
}

# Function to plot heatmap ordered by p(clusterK)
plotOrderedHeat <- function (cluster = "k1") {
  ids_order <- as.character(clust_res$eid[orderIDs(cluster)])
  # Arrange B-coef matrix
  plot_mat <- as.matrix(mn_mat[ids_order, ])
  # For plotting, scale each column (RINT) to see differences better
  plot_mat <- apply(plot_mat, 2, FUN = RINTx)
  
  png(paste0(plotdir, "heatmaps_ordered_", cluster, ".png"),
      width = 21, height = 29.7, units = "cm", res = 300)
  pheatmap(plot_mat, 
           cluster_rows = F, cluster_cols = F,
           scale = "none", 
           treeheight_col = 0, treeheight_row = 0,
           show_rownames = F, show_colnames = F,
           color = colorRampPalette(c("navy", "white", "red"))(50),
           legend = T, fontsize = 10,
           main = paste0("B-coefs in ", PHENO, "_", SEX_STRATA))
  dev.off()
}

# Probability bar plots to go alongside the histogram
plotOrderedBars <- function (cluster = "k1") {
  # Rename IDs and flip the factor levels to retain correct order for plot
  dat <- clust_res[orderIDs(cluster), ]
  dat$eid <- factor(as.character(1:nrow(dat)),
                    levels = as.character(nrow(dat):1))
  
  resplot <- ggplot(dat, aes(x = !!as.symbol(cluster),
                             y = eid)) +
    geom_col(width = 1, 
             color = custom_four_diverge[cluster],
             fill = custom_four_diverge[cluster]) +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank())
  png(paste0(plotdir, "barprobs_", cluster, ".png"),
      width = 21/2, height = 29.7, units = "cm", res = 300)
  print(resplot)
  dev.off()
}

lapply(paste0("k", 1:K_chosen), function (clustk) {
  plotOrderedHeat(clustk)
  plotOrderedBars(clustk)
})

## Population trajectories ----

# On the scale of adjusted values and time from first measurement
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
    scale_color_manual(values = custom_four_diverge) +
    scale_fill_manual(values = custom_four_diverge) +
    theme(legend.position = "none") +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE))
  return (popn_traj_plot)
}

# Apply to adj- and baselined- data

to_plot_adj_baselined <- lapply(paste0("k", 1:K_chosen), function (clustk) {
  summ_clust <- getSummDat(cluster = clustk, 
                           prob_assign = 0.75,
                           valtype = "value_fulladj", 
                           timetype = "t_diff") %>%
    mutate(clust = clustk)
  return (summ_clust)
})
to_plot_adj_baselined <- bind_rows(to_plot_adj_baselined) %>%
  mutate(clust = factor(clust, 
                        levels = paste0("k", 1:K_chosen)))
png(paste0(plotdir, "popn_trajectories_adj_dat.png"))
print(plotTrajectories(to_plot_adj_baselined))
dev.off()

# Apply to raw data

to_plot_raw <- lapply(paste0("k", 1:K_chosen), function (clustk) {
  summ_clust <- getSummDat(cluster = clustk, 
                           prob_assign = 0.75,
                           valtype = "value", 
                           timetype = "age_t1") %>%
    mutate(clust = clustk)
  return (summ_clust)
})
to_plot_raw <- bind_rows(to_plot_raw) %>%
  mutate(clust = factor(clust, 
                        levels = paste0("k", 1:K_chosen)))
png(paste0(plotdir, "popn_trajectories_raw.png"))
print(plotTrajectories(to_plot_raw) +
        scale_x_continuous(limits = c(30, 70)))
dev.off()
