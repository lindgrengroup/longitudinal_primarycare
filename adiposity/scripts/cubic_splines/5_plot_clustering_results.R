# Author: Samvida S. Venkatesh
# Date: 14/12/21

library(lme4)
library(splines)
library(tidyverse)
library(zoo)
library(RColorBrewer)
theme_set(theme_bw())

set.seed(141221)

# Read files ----

PHENOTYPES <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/pheno_names.txt")$V1
SEX_STRATA <- c("F", "M", "sex_comb")

slope_models <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines_",
                 p, ".rds"))
})
names(slope_models) <- PHENOTYPES

cluster_membership <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SEX_STRATA, function (sx) {
    read.table(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_spline_clusters_",
                      p, "_", sx, ".txt"),
               sep = "\t", header = T, stringsAsFactors = F)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(cluster_membership) <- PHENOTYPES

dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/indiv_qcd_data.rds")[PHENOTYPES]
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/covariates.rds")[PHENOTYPES]

PCs <- paste0("PC", 1:21)
COVARS <- c("baseline_age", "age_sq", "data_provider")

# Wrangle data ----

# Add in covariates to raw GP data
model_dat <- lapply(PHENOTYPES, function (p) {
  res <- merge(dat[[p]], covars[[p]], by = "eid")
  return (res)
})
names(model_dat) <- PHENOTYPES

# Count number of individuals in each cluster to make sampling quicker later
cluster_membership <- lapply(cluster_membership, function (list_df) {
  lapply(list_df, function (df) {
    res <- df %>% group_by(k) %>% mutate(n_cluster = n())
    return (res)
  })
})

# Plot raw data and model predictions for random individuals ----
## Get random IDs within each cluster ----

get_rand_eids <- function (p, sx, n_each = 5) {
  # Return df of id and cluster number
  sampled_ids <- cluster_membership[[p]][[sx]] %>% 
    group_by(k) %>%
    sample_n(ifelse(n_cluster < n_each, n_cluster, n_each)) %>%
    select(eid, k)
  return (sampled_ids)
}

## Function to create predicted data for set of ids ----

create_prediction_df <- function (p, sx, ids) {
  sub_dat <- subset(model_dat[[p]], 
                    model_dat[[p]]$eid %in% ids)
  # Calculate maximum time-point to predict to for each eid
  sub_dat <- sub_dat %>% group_by(eid) %>% 
    mutate(min_age = floor(min(age_event)),
           max_age = ceiling(max(age_event)))
  
  # Create new data to predict from
  new_data <- sub_dat %>% select(all_of(c("eid", COVARS, PCs,
                                          "min_age", "max_age", "sex"))) %>% 
    distinct(eid, data_provider, .keep_all = T) 
  # Timepoints to extend to (create prediction every 3 months) 
  per_id <- lapply(1:nrow(new_data), FUN = function (i) { 
    seq(new_data$min_age[i], new_data$max_age[i], by = 0.25) 
  })
  length_per_id <- unlist(lapply(per_id, function (x) length (x)))
  t <- unlist(per_id)
  tmp <- data.frame(eid = rep(new_data$eid, length_per_id),
                    data_provider = rep(new_data$data_provider, length_per_id),
                    age_event = t)
  new_data <- merge(tmp, new_data, by = c("eid", "data_provider"))
  
  # Predict new values
  fitted_results <- 
    as.data.frame(predict(slope_models[[p]][[sx]],
                          newdata = new_data))
  colnames(fitted_results) <- "fit"
  pred_df <- bind_cols(new_data, fitted_results)
  
  # At each time-point, average across data providers for each individual
  pred_df <- pred_df %>% group_by(eid, age_event) %>%
    summarise(fit = mean(fit))
  
  return (pred_df)
}

## Function to create plots ----

plot_predictions <- function (p, sx, id_df) {
  raw_dat <- model_dat[[p]] %>% filter(eid %in% id_df$eid)
  
  plot_dat <- create_prediction_df(p, sx, id_df$eid)
  plot_dat$cluster <- id_df$k[match(plot_dat$eid,
                                    id_df$eid)]
  plot_dat$cluster <- as.factor(as.character(plot_dat$cluster))
  
  # Order the ids by cluster number for plot
  id_levels <- id_df$eid[order(id_df$k)]
  plot_dat$eid_f <- factor(as.character(plot_dat$eid), levels = id_levels)
  raw_dat$eid_f <- factor(as.character(raw_dat$eid), levels = id_levels)
  
  res <- ggplot(plot_dat, aes(x = age_event)) +
    facet_wrap(.~eid_f, ncol = 5, scales = "free") +
    geom_point(data = raw_dat, aes(y = value),
               colour = "black") +
    geom_line(aes(y = fit, colour = cluster)) +
    scale_color_brewer(palette = "Set1") +
    labs(x = "Age (years)",
         y = p,
         title = paste0("Randomly selected ", sx, ", phenotype: ", p))
  return (res)
}

# Plot model predictions for hundreds of samples ----

plot_samples_modelled <- function (p, sx) {
  # Get up to 100 random IDs within each cluster
  ids_to_plot <- get_rand_eids(p, sx, 100)
  # Within each cluster, choose 5 random ids to highlight
  highlight_ids <- ids_to_plot %>% group_by(k) %>% sample_n(5)
  highlight_ids <- unique(highlight_ids$eid)
  
  # Create fitted data
  plot_dat <- create_prediction_df(p, sx, ids_to_plot$eid)
  plot_dat$cluster <- ids_to_plot$k[match(plot_dat$eid,
                                          ids_to_plot$eid)]
  plot_dat$cluster <- as.factor(as.character(plot_dat$cluster))
  plot_dat$highlight <- plot_dat$eid %in% highlight_ids
  
  # Plot
  res_plot <- ggplot(plot_dat %>% filter(!highlight), 
                     aes(x = age_event, y = fit, group = eid,
                         colour = cluster)) +
    facet_wrap(~cluster) +
    geom_line(alpha = 0.1) +
    geom_line(data = plot_dat %>% filter(highlight),
              alpha = 1) +
    scale_color_brewer(palette = "Set1") +
    labs(x = "Age (years)",
         y = p,
         title = paste0("Sample of modelled trajectories in each cluster of phenotype: ",
                        p, " strata: ", sx))
  return (res_plot)
}

# Plot mean model prediction in each cluster ----

plot_all_modelled <- function (p, sx) {
  
  # Get list of individuals in strata
  cluster_info <- cluster_membership[[p]][[sx]]
  clusts <- unique(cluster_info$k)
  # Get predictions in each cluster
  per_clust_pred <- lapply(clusts, function (k) {
    clust_ids <- cluster_info$eid[cluster_info$k == k]
    pred_df <- create_prediction_df(p, sx, clust_ids)
    pred_df$cluster <- k
    return (pred_df)
  })
  # Bind results across all clusters
  pred_dat <- bind_rows(per_clust_pred)
  
  # At each time-point, average across all fitted data to get mean and s.e.
  # of mean (95% CI) within each cluster
  plot_dat <- pred_dat %>% group_by(cluster, age_event) %>%
    summarise(mean_fit = mean(fit),
              sd_fit = sd(fit), 
              n = n()) %>%
    mutate(lci = mean_fit - 1.96*sd_fit,
           uci = mean_fit + 1.96*sd_fit,
           lci_mean = mean_fit - 1.96*(sd_fit/sqrt(n)),
           uci_mean = mean_fit + 1.96*(sd_fit/sqrt(n))) 
  # Get rolling average of these stats across 2 years 
  # because the plots are too choppy otherwise
  fn_roll <- function (x) { rollapply(x, 8, mean, fill = NA) }
  plot_dat <- plot_dat %>% ungroup() %>% 
    group_by(cluster) %>% arrange(age_event, .by_group = T) %>%
    mutate(across(c(mean_fit, lci, uci, lci_mean, uci_mean),
                  fn_roll))
  
  # Plot
  plot_dat$cluster <- as.factor(as.character(plot_dat$cluster))
  
  facet_plot <- ggplot(plot_dat,
                       aes(x = age_event, y = mean_fit, 
                           fill = cluster, 
                           colour = cluster)) +
    facet_wrap(~cluster) +
    geom_line() +
    geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.1) +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.3) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Age (years)", 
         y = paste0("Mean and 95% C.I. of ", p), 
         title = paste0("Modelled trajectories in each cluster of phenotype: ",
                        p, " strata: ", sx))
  
  all_plot <- ggplot(plot_dat,
                     aes(x = age_event, y = mean_fit, 
                         fill = cluster, 
                         colour = cluster)) +
    geom_line() +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.1) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Age (years)", 
         y = paste0("Mean and 95% C.I. of mean of ", p), 
         title = paste0("Modelled trajectories in each cluster of phenotype: ",
                        p, " strata: ", sx))
  
  return (list(facet_plot, all_plot))
}

# Plot population-level trajectories ----

plot_mean_traj_by_cluster <- function (p, sx) {
  raw_dat_with_k <- model_dat[[p]] %>%
    filter(eid %in% cluster_membership[[p]][[sx]]$eid) 
  raw_dat_with_k$cluster <- 
    cluster_membership[[p]][[sx]]$k[match(raw_dat_with_k$eid,
                                          cluster_membership[[p]][[sx]]$eid)]
  # Cut ages into 1 year bins
  cut_points <- seq(20, 80, by = 1)
  cut_labels <- seq(20, 79, by = 1)
  raw_dat_with_k$age_bin <- cut(raw_dat_with_k$age_event, 
                                breaks = cut_points, 
                                include.lowest = T,
                                labels = cut_labels)
  summ_dat <- raw_dat_with_k %>% 
    group_by(cluster, age_bin) %>% 
    summarise(mean_value = mean(value),
              sd_value = sd(value), 
              n = n()) %>%
    mutate(lci = mean_value - 1.96*sd_value,
           uci = mean_value + 1.96*sd_value,
           lci_mean = mean_value - 1.96*(sd_value/sqrt(n)),
           uci_mean = mean_value + 1.96*(sd_value/sqrt(n))) 
  
  # Get rolling average of these stats across 5 years 
  # because the plots are too choppy otherwise
  fn_roll <- function (x) { rollapply(x, 5, mean, fill = NA) }
  summ_dat <- summ_dat %>% ungroup() %>% 
    group_by(cluster) %>% mutate(age_bin = as.numeric(age_bin)) %>%
    arrange(age_bin, .by_group = T) %>%
    mutate(across(c(mean_value, lci, uci, lci_mean, uci_mean),
                  fn_roll))
  
  summ_dat$cluster <- as.factor(as.character(summ_dat$cluster))
  
  facet_plot <- ggplot(summ_dat, aes(x = age_bin, y = mean_value, 
                                   color = cluster, fill = cluster)) +
    facet_wrap(~cluster) +
    geom_line() +
    geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.1) +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.3) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Age (years)", 
         y = paste0("Mean and 95% C.I. of ", p), 
         title = paste0("Observed adiposity in each cluster of phenotype: ",
                        p, " strata: ", sx))
  
  all_plot <- ggplot(summ_dat,
                     aes(x = age_bin, y = mean_value, 
                         color = cluster, fill = cluster)) +
    geom_line() +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.1) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Age (years)", 
         y = paste0("Mean and 95% C.I. of mean of ", p), 
         title = paste0("Modelled trajectories in each cluster of phenotype: ",
                        p, " strata: ", sx))
  
  return (list(facet_plot, all_plot))
  
}

# Plot cluster identity vs covariates ----

sex_col_palette <- c("#F8766D", "#00BFC4", "#C77CFF")
names(sex_col_palette) <- c("F", "M", "sex_comb")

plot_cluster_vs_covs <- function (p, sx) {
  # Get covariates data
  covars_with_k <- covars[[p]] %>%
    filter(eid %in% cluster_membership[[p]][[sx]]$eid) 
  covars_with_k$cluster <- 
    cluster_membership[[p]][[sx]]$k[match(covars_with_k$eid,
                                          cluster_membership[[p]][[sx]]$eid)]
  covars_with_k$cluster <- as.factor(as.character(covars_with_k$cluster))
  
  # Number of men/women (individuals) in each cluster
  summ_sex <- covars_with_k %>% count(cluster, sex)
  sex_plot <- ggplot(summ_sex, 
                     aes(x = cluster, y = n, fill = sex)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(values = sex_col_palette) +
    labs(x = "Cluster identity", y = "Number of individuals", 
         title = paste0("Sex distribution in clusters of: ",
                        p, " strata: ", sx)) 
  
  # Baseline age stratified by cluster identity
  bl_age_plot <- ggplot(covars_with_k, 
                        aes(x = cluster, y = baseline_age)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Set1") + 
    labs(x = "Cluster identity", y = "Baseline age (years)", 
         title = paste0("Baseline age distribution in clusters of: ",
                        p, " strata: ", sx)) +
    theme(legend.position = "none")
  
  # Baseline trait value
  bl_trait_plot <- ggplot(covars_with_k, 
                          aes(x = cluster, y = baseline_trait)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Set1") + 
    labs(x = "Cluster identity", y = "Baseline phenotype", 
         title = paste0("Baseline phenotype distribution in clusters of: ",
                        p, " strata: ", sx)) +
    theme(legend.position = "none")
  
  # Number of follow-up measures
  max_y <- min(50, max(covars_with_k$FU_n))
  fu_n_plot <- ggplot(covars_with_k, 
                      aes(x = cluster, y = FU_n)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Set1") + 
    scale_y_continuous(limits = c(0, max_y)) +
    labs(x = "Cluster identity", y = "# follow-up measures", 
         title = paste0("Distribution of # follow-ups in clusters of: ",
                        p, " strata: ", sx)) +
    theme(legend.position = "none")
  
  # Number of follow-up years
  fu_yrs_plot <- ggplot(covars_with_k, 
                        aes(x = cluster, y = FUyrs)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Set1") + 
    labs(x = "Cluster identity", y = "# years follow-up", 
         title = paste0("Distribution of # years of follow-up in clusters of: ",
                        p, " strata: ", sx)) +
    theme(legend.position = "none")
  
  return (list(sex_plot, 
               bl_age_plot, bl_trait_plot, 
               fu_n_plot, fu_yrs_plot))
}

# Apply all plotting functions ----

plot_list <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    
    # Modelled and observed trajectories of random ids
    rand_plots <- plot_predictions(p, sx, get_rand_eids(p, sx, 5))
    # Modelled trajectories of random ids by cluster
    cluster_sample_plots <- plot_samples_modelled(p, sx)
    # Modelled trajectories (all in cluster)
    modelled_all_plots <- plot_all_modelled(p, sx)
    # Population-level trajectories
    popn_plots <- plot_mean_traj_by_cluster(p, sx)
    # Cluster membership vs covariates 
    cov_plots <- plot_cluster_vs_covs(p, sx)
    
    pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/cubic_splines/clustering/cluster_trajectories_",
               p, "_", sx, ".pdf"))
    print(rand_plots)
    print(cluster_sample_plots)
    print(modelled_all_plots)
    print(popn_plots)
    print(cov_plots)
    dev.off()
    
    return ()
  })
  return ()
})

