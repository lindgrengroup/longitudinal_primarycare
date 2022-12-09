# Author: Samvida S. Venkatesh
# Date: 24/05/2022

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

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")
K_chosen <- 4

plotdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/highdim_splines_clustering/plots/"
resdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity//ukb_no_gp/highdim_splines_clustering/"

# Load data ----

# Soft clustering probabilities
clust_res <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    res <- read.table(paste0(resdir, p, "_", sx, "/soft_clustering_probs_", p, "_", sx, 
                      ".txt"),
               sep = "\t", header = T, stringsAsFactors = F)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(clust_res) <- PHENOTYPES

# Covariates
general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

# Original data
orig_popn_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/highdim_splines_clustering/data/dat_to_model.rds")

# Plot soft clustering results ----

CAT_COVARS <- "sex"
QUANT_COVARS <- c("baseline_age", "baseline_trait",
                  "FUyrs", "FU_n",
                  "year_of_birth", "age_at_death")

# Create covars
covars <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    res <- orig_popn_dat[[p]][[sx]] %>%
      group_by(eid) %>% arrange(age_t1, .by_group = T) %>%
      mutate(baseline_age = first(age_t1),
             baseline_trait = first(value_t1),
             FUyrs = (t_diff - 1)/365,
             FU_n = n()) %>%
      # Only keep last occurrence as that gives us # years and length of follow-up
      slice(tail(row_number(), 1)) %>%
      select(any_of(c("eid", CAT_COVARS, QUANT_COVARS)))
    
    res <- left_join(res, general_covars, by = "eid")
    return (res)
  })
  names(per_sex) <- SEX_STRATA
  return (per_sex)
})
names(covars) <- PHENOTYPES

## Plots for max probability per individual -----

# Histogram
lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    dir.create(paste0(plotdir, p, "_", sx))
    for_hist <- clust_res[[p]][[sx]] %>%
      mutate(max_prob = dplyr::select(., starts_with("k")) %>% 
               do.call(pmax, .))
    
    hist_plot <- ggplot(for_hist, aes(x = max_prob)) +
      geom_histogram() +
      expand_limits(x = 0)
    
    png(paste0(plotdir, p, "_", sx, "/max_probabilities_hist.png"))
    print(hist_plot)
    dev.off()
  })
})

# vs covariates
lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    for_hist <- clust_res[[p]][[sx]] %>%
      mutate(max_prob = dplyr::select(., starts_with("k")) %>% 
               do.call(pmax, .)) %>%
      dplyr::select(all_of(c("eid", "max_prob"))) %>%
      # Cut max prob into 0.1 unit bins (<0.5, <0.6, etc.)
      mutate(max_prob = as.factor(as.character(paste0("<", 
                                                      plyr::round_any(max_prob, 0.1, f = floor)))),
             eid = as.character(eid))
    
    covar_hists <- left_join(for_hist, covars[[p]][[sx]], by = "eid")
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
    png(paste0(plotdir, p, "_", sx, "/max_probabilities_vs_covars.png"))
    print(arranged_plots)
    dev.off()
  })
})

## Plots for probability of one cluster vs another -----

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    pairs_plot <- ggpairs(clust_res[[p]][[sx]], columns = 2:(K_chosen+1))
    png(paste0(plotdir, p, "_", sx, "/clust_probs_pairs_plot.png"))
    print(pairs_plot)
    dev.off()
  })
})

## Population trajectories ----

# On the scale of adjusted values and time from first measurement
getSummDat <- function (clust_dat, popn_dat,
                        cluster = "k1", prob_assign = 0.75,
                        valtype = "value", timetype = "age_t1") {
  keep_ids <- as.character(clust_dat$eid[which(clust_dat[, cluster] >= prob_assign)])
  
  if (timetype == "age_t1") {
    round_yrs = 0.5
    roll_yrs = 2
  } else if (timetype == "t_diff") {
    round_yrs = 500
    roll_yrs = 2000
  }
  
  dat_summ <- popn_dat %>%
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

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    to_plot_adj_baselined <- lapply(paste0("k", 1:K_chosen), function (clustk) {
      summ_clust <- getSummDat(clust_dat = clust_res[[p]][[sx]],
                               popn_dat = orig_popn_dat[[p]][[sx]],
                               cluster = clustk, 
                               prob_assign = 0.5,
                               valtype = "value_fulladj", 
                               timetype = "t_diff") %>%
        mutate(clust = clustk)
      return (summ_clust)
    })
    to_plot_adj_baselined <- bind_rows(to_plot_adj_baselined) %>%
      mutate(clust = factor(clust, 
                            levels = paste0("k", 1:K_chosen)))
    png(paste0(plotdir, p, "_", sx, "/popn_trajectories_adj_dat.png"))
    print(plotTrajectories(to_plot_adj_baselined) +
            scale_x_continuous(limits = c(0, 2500)))
    dev.off()
  })
})

# Apply to raw data

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    to_plot_raw <- lapply(paste0("k", 1:K_chosen), function (clustk) {
      summ_clust <- getSummDat(clust_dat = clust_res[[p]][[sx]],
                               popn_dat = orig_popn_dat[[p]][[sx]],
                               cluster = clustk, 
                               prob_assign = 0.5,
                               valtype = "value", 
                               timetype = "age_t1") %>%
        mutate(clust = clustk)
      return (summ_clust)
    })
    to_plot_raw <- bind_rows(to_plot_raw) %>%
      mutate(clust = factor(clust, 
                            levels = paste0("k", 1:K_chosen)))
    png(paste0(plotdir, p, "_", sx, "/popn_trajectories_raw.png"))
    print(plotTrajectories(to_plot_raw) +
            scale_x_continuous(limits = c(40, 70)))
    dev.off()
  })
})

