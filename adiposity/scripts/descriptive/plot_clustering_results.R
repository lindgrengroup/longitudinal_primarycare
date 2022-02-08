# Author: Samvida S. Venkatesh
# Date: 14/12/21

library(lme4)
library(splines)
library(tidyverse)
library(zoo)
library(pheatmap)
library(RColorBrewer)
theme_set(theme_bw())

set.seed(141221)

# Read files ----

args <- commandArgs(trailingOnly = T)
STRATA <- args[1]

p <- gsub("_.*", "", STRATA)
sx <- gsub(paste0(p, "_"), "", STRATA)

slope_models <- readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/age_adj_",
                               p, ".rds"))[[sx]]

blups <- readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/age_adj_blups_",
                        p, ".rds"))[[sx]]

cluster_membership <- read.table(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/age_adj_clusters_",
                                        p, "_", sx, ".txt"),
                                 sep = "\t", header = T, stringsAsFactors = F)

dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/indiv_qcd_data.rds")[[p]]
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/covariates.rds")[[p]]

NFE <- ifelse(p %in% c("BMI", "Weight"), 3, 1) 
NRE <- ifelse(p %in% c("BMI", "Weight"), 3, 1) 

# Wrangle data ----

# Add in covariates to raw GP data
model_dat <- merge(dat, covars, by = "eid")
if (sx != "sex_comb") model_dat <- model_dat %>% filter(sex == sx)
model_dat <- model_dat %>% mutate(t = age_event - baseline_age)

# Count number of individuals in each cluster to make sampling quicker later
cluster_membership <- cluster_membership %>% 
  group_by(k) %>% 
  mutate(n_cluster = n())

# Plot raw data and model predictions for random individuals ----
## Get random IDs within each cluster ----

get_rand_eids <- function (n_each = 5) {
  # Return df of id and cluster number
  sampled_ids <- cluster_membership %>% 
    group_by(k) %>%
    sample_n(ifelse(n_cluster < n_each, n_cluster, n_each)) %>%
    select(eid, k)
  return (sampled_ids)
}

## Function to create predicted data for set of ids ----

create_prediction_df <- function (model_name, ids) {
  sub_dat <- subset(model_dat, 
                    model_dat$eid %in% ids)
  # Calculate maximum time-point to predict to for each eid
  sub_dat <- sub_dat %>% group_by(eid) %>% 
    mutate(max_t = ceiling(max(t)))
  
  # Create new data to predict from
  new_data <- sub_dat %>% 
    select(all_of(c("eid", "baseline_age", "age_sq", "max_t", "sex"))) %>% 
    distinct(eid, .keep_all = T) 
  # Timepoints to extend to
  ts <- lapply(1:nrow(new_data), FUN = function (i) { 
    seq(0, new_data$max_t[i], by = 0.25) })
  length_ts <- unlist(lapply(ts, function (x) length(x)))
  t <- unlist(ts)
  tmp <- data.frame(eid = rep(new_data$eid, length_ts),
                    t = t)
  new_data <- merge(tmp, new_data, by = "eid")
  
  # Predict new values
  fitted_results <- 
    as.data.frame(predict(model_name, newdata = new_data))
  colnames(fitted_results) <- "fit"
  pred_df <- bind_cols(new_data, fitted_results)
  
  # Add column for age at event
  pred_df <- pred_df %>% 
    mutate(age_event = t + baseline_age) %>%
    select(eid, t, age_event, fit)
  
  return (pred_df)
}

## Function to create plots ----

plot_predictions <- function (id_df) {
  raw_dat <- model_dat %>% filter(eid %in% id_df$eid)
  
  plot_dat <- create_prediction_df(slope_models, id_df$eid)
  plot_dat$cluster <- id_df$k[match(plot_dat$eid,
                                    id_df$eid)]
  plot_dat$cluster <- as.factor(as.character(plot_dat$cluster))
  
  # Order the ids by cluster number for plot
  id_levels <- id_df$eid[order(id_df$k)]
  plot_dat$eid_f <- factor(as.character(plot_dat$eid), levels = id_levels)
  raw_dat$eid_f <- factor(as.character(raw_dat$eid), levels = id_levels)
  
  res <- ggplot(plot_dat, aes(x = age_event)) +
    facet_wrap(.~eid_f, ncol = 5, scales = "free_x") +
    geom_point(data = raw_dat, aes(y = value),
               colour = "black") +
    geom_line(aes(y = fit, colour = cluster)) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Age (years)",
         y = p,
         title = paste0("Randomly selected ", sx, ", phenotype: ", p))
  return (res)
}

# Plot model predictions for hundreds of samples ----

plot_samples_modelled <- function () {
  # Get up to 100 random IDs within each cluster
  ids_to_plot <- get_rand_eids(100)
  # Within each cluster, choose 5 random ids to highlight
  highlight_ids <- ids_to_plot %>% group_by(k) %>% sample_n(5)
  highlight_ids <- unique(highlight_ids$eid)
  
  # Create fitted data
  plot_dat <- create_prediction_df(slope_models, ids_to_plot$eid)
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
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Age (years)",
         y = p,
         title = paste0("Sample of modelled trajectories in each cluster of phenotype: ",
                        p, " strata: ", sx))
  return (res_plot)
}

# Re-fit models and plot fixed effect in each cluster ----

plot_refit_clusters <- function () {
  # Create modelled data split by cluster
  to_model <- model_dat 
  
  to_model$k <- 
    cluster_membership$k[match(to_model$eid,
                               cluster_membership$eid)]
  to_model <- split(to_model, f = to_model$k, drop = F)
  
  fitted_dat_within_cluster <- lapply(1:length(to_model), function (ki) {
    df <- to_model[[ki]]
    mod_covars <- c("baseline_age", "age_sq")
    if (sx == "sex_comb") mod_covars <- c(mod_covars, "sex")
    mod_formula <- formula(paste0("value ~ ", 
                                  paste0(mod_covars, collapse = " + "), 
                                  "+ ns(t, df = ", NFE, ") + (ns(t, df = ", NRE, ") | eid)"))
    mod_res <- lmer(mod_formula, data = df, REML = F)
    # Fixed effect prediction 
    newdat <- create_prediction_df(mod_res, df$eid)
    newdat_sumstats <- newdat %>% group_by(t) %>% 
      summarise(mean_fit = mean(fit),
                sd_fit = sd(fit), 
                n = n()) %>%
      mutate(lci = mean_fit - 1.96*sd_fit,
             uci = mean_fit + 1.96*sd_fit,
             lci_mean = mean_fit - 1.96*(sd_fit/sqrt(n)),
             uci_mean = mean_fit + 1.96*(sd_fit/sqrt(n)),
             cluster = ki) 
    return (newdat_sumstats)
  })
  plot_dat <- bind_rows(fitted_dat_within_cluster)
  
  # Get rolling average of these stats across 2 years 
  # because the plots are too choppy otherwise
  fn_roll <- function (x) { rollapply(x, 8, mean, fill = NA) }
  plot_dat <- plot_dat %>% ungroup() %>% 
    group_by(cluster) %>% arrange(t, .by_group = T) %>%
    mutate(across(c(mean_fit, lci, uci, lci_mean, uci_mean),
                  fn_roll))
  
  # Plot
  plot_dat$cluster <- as.factor(as.character(plot_dat$cluster))
  
  facet_plot <- ggplot(plot_dat,
                       aes(x = t, y = mean_fit, 
                           fill = cluster, 
                           colour = cluster)) +
    facet_wrap(~cluster) +
    geom_line() +
    geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.1) +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.3) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "Time from baseline measurement (years)", 
         y = paste0("Mean and 95% C.I. of ", p), 
         title = paste0("Modelled trajectories in each cluster of phenotype: ",
                        p, " strata: ", sx))
  
  all_plot <- ggplot(plot_dat,
                     aes(x = t, y = mean_fit, 
                         fill = cluster, 
                         colour = cluster)) +
    geom_line() +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.1) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "Time from baseline measurement (years)", 
         y = paste0("Mean and 95% C.I. of mean of ", p), 
         title = paste0("Modelled trajectories in each cluster of phenotype: ",
                        p, " strata: ", sx))
  
  return (list(facet_plot, all_plot))
  
}

# # Plot mean model prediction in each cluster ----
# 
# plot_all_modelled <- function (p, sx) {
#   
#   # Get list of individuals in strata
#   cluster_info <- cluster_membership[[p]][[sx]]
#   clusts <- unique(cluster_info$k)
#   # Get predictions in each cluster
#   per_clust_pred <- lapply(clusts, function (k) {
#     clust_ids <- cluster_info$eid[cluster_info$k == k]
#     pred_df <- create_prediction_df(slope_models[[p]][[sx]], p, sx, clust_ids)
#     pred_df$cluster <- k
#     return (pred_df)
#   })
#   # Bind results across all clusters
#   pred_dat <- bind_rows(per_clust_pred)
#   
#   # At each time-point, average across all fitted data to get mean and s.e.
#   # of mean (95% CI) within each cluster
#   plot_dat <- pred_dat %>% ungroup() %>%
#     group_by(cluster, t) %>% summarise(mean_fit = mean(fit),
#                                        sd_fit = sd(fit), 
#                                        n = n()) %>%
#     mutate(lci = mean_fit - 1.96*sd_fit,
#            uci = mean_fit + 1.96*sd_fit,
#            lci_mean = mean_fit - 1.96*(sd_fit/sqrt(n)),
#            uci_mean = mean_fit + 1.96*(sd_fit/sqrt(n))) 
#   # Get rolling average of these stats across 2 years 
#   # because the plots are too choppy otherwise
#   fn_roll <- function (x) { rollapply(x, 8, mean, fill = NA) }
#   plot_dat <- plot_dat %>% ungroup() %>% 
#     group_by(cluster) %>% arrange(t, .by_group = T) %>%
#     mutate(across(c(mean_fit, lci, uci, lci_mean, uci_mean),
#                   fn_roll))
#   
#   # Plot
#   plot_dat$cluster <- as.factor(as.character(plot_dat$cluster))
#   
#   facet_plot <- ggplot(plot_dat,
#                        aes(x = t, y = mean_fit, 
#                            fill = cluster, 
#                            colour = cluster)) +
#     facet_wrap(~cluster) +
#     geom_line() +
#     geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.1) +
#     geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.3) +
#     scale_color_brewer(palette = "Dark2") +
#     scale_fill_brewer(palette = "Dark2") +
#     labs(x = "Time from baseline measurement (years)", 
#          y = paste0("Mean and 95% C.I. of ", p), 
#          title = paste0("Modelled trajectories in each cluster of phenotype: ",
#                         p, " strata: ", sx))
#   
#   all_plot <- ggplot(plot_dat,
#                      aes(x = t, y = mean_fit, 
#                          fill = cluster, 
#                          colour = cluster)) +
#     geom_line() +
#     geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.1) +
#     scale_color_brewer(palette = "Dark2") +
#     scale_fill_brewer(palette = "Dark2") +
#     labs(x = "Time from baseline measurement (years)", 
#          y = paste0("Mean and 95% C.I. of mean of ", p), 
#          title = paste0("Modelled trajectories in each cluster of phenotype: ",
#                         p, " strata: ", sx))
#   
#   return (list(facet_plot, all_plot))
# }

# Plot population-level observed trajectories ----

plot_mean_traj_by_cluster <- function () {
  raw_dat_with_k <- model_dat %>%
    filter(eid %in% cluster_membership$eid) 
  raw_dat_with_k$cluster <- 
    cluster_membership$k[match(raw_dat_with_k$eid,
                               cluster_membership$eid)]
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
    group_by(cluster) %>% 
    mutate(age_bin = as.numeric(as.character(age_bin))) %>%
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
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "Age (years)", 
         y = paste0("Mean and 95% C.I. of ", p), 
         title = paste0("Observed adiposity in each cluster of phenotype: ",
                        p, " strata: ", sx))
  
  all_plot <- ggplot(summ_dat,
                     aes(x = age_bin, y = mean_value, 
                         color = cluster, fill = cluster)) +
    geom_line() +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.1) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "Age (years)", 
         y = paste0("Mean and 95% C.I. of mean of ", p), 
         title = paste0("Observed trajectories in each cluster of phenotype: ",
                        p, " strata: ", sx))
  
  return (list(facet_plot, all_plot))
  
}

# Plot cluster identity vs covariates ----

sex_col_palette <- c("#F8766D", "#00BFC4", "#C77CFF")
names(sex_col_palette) <- c("F", "M", "sex_comb")

plot_cluster_vs_covs <- function () {
  # Get covariates data
  covars_with_k <- covars %>%
    filter(eid %in% cluster_membership$eid) 
  covars_with_k$cluster <- 
    cluster_membership$k[match(covars_with_k$eid,
                               cluster_membership$eid)]
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
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "Cluster identity", y = "Baseline age (years)", 
         title = paste0("Baseline age distribution in clusters of: ",
                        p, " strata: ", sx)) +
    theme(legend.position = "none")
  
  # Baseline trait value
  bl_trait_plot <- ggplot(covars_with_k, 
                          aes(x = cluster, y = baseline_trait)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Dark2") + 
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
    scale_fill_brewer(palette = "Dark2") + 
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
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "Cluster identity", y = "# years follow-up", 
         title = paste0("Distribution of # years of follow-up in clusters of: ",
                        p, " strata: ", sx)) +
    theme(legend.position = "none")
  
  return (list(sex_plot, 
               bl_age_plot, bl_trait_plot, 
               fu_n_plot, fu_yrs_plot))
}

# Apply all plotting functions ----

# Modelled and observed trajectories of random ids
rand_plots <- plot_predictions(get_rand_eids(5))

# Modelled trajectories of random ids by cluster
cluster_sample_plots <- plot_samples_modelled()

# Refit trajectories within each cluster
refit_plots <- plot_refit_clusters()

# Modelled trajectories (all in cluster)
# modelled_all_plots <- plot_all_modelled(p, sx)

# Population-level observed trajectories
popn_plots <- plot_mean_traj_by_cluster()

# Cluster membership vs covariates 
cov_plots <- plot_cluster_vs_covs()

pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/cubic_splines/clustering/cluster_trajectories_",
           p, "_", sx, ".pdf"))
print(rand_plots)
print(cluster_sample_plots)
print(refit_plots)
# print(modelled_all_plots)
print(popn_plots)
print(cov_plots)
dev.off()

# Plot heatmaps of spline coefficients in the clusters ----

# RINT coefficient for plotting
rint_x <- function (x) {
  return (qnorm((rank(x) - 0.5) / sum(!is.na(x))))
}

# Row annotations (cluster membership)
k_info <- cluster_membership
# Order BLUPs by cluster membership
k_info <- k_info[order(k_info$k), ]
row_annots <- data.frame(cluster = as.character(k_info$k))
rownames(row_annots) <- as.character(k_info$eid)

for_plot <- blups[as.character(k_info$eid), ]
for_plot <- as.matrix(for_plot)
int_col <- grep("Intercept", colnames(for_plot))
# Replace parentheses text with blanks but retain (Intercept)
colnames(for_plot) <- gsub("\\s*\\([^\\)]+\\)", "", 
                           colnames(for_plot))
colnames(for_plot)[int_col] <- "Intercept"

# For plotting, scale each column (RINT) to see differences better
for_plot <- apply(for_plot, 2, FUN = rint_x)

pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/cubic_splines/clustering/cluster_heatmaps_",
           p, "_", sx, ".pdf"))
pheatmap(for_plot, 
         cluster_rows = F, cluster_cols = F,
         scale = "none",   
         annotation_row = row_annots, 
         treeheight_col = 0, treeheight_row = 50,
         show_rownames = F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         legend = T, fontsize = 10,
         main = paste0("Cubic spline BLUPs in ", p, "_", sx))
dev.off()

