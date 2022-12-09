# Author: Samvida S. Venkatesh
# Date: 14/12/21

library(lme4)
library(splines)
library(tidyverse)
library(zoo)
library(pheatmap)
library(ggpubr)
library(RColorBrewer)
theme_set(theme_bw())

set.seed(141221)

# Read files ----

args <- commandArgs(trailingOnly = T)
STRATA <- args[1]

p <- gsub("_.*", "", STRATA)
sx <- gsub(paste0(p, "_"), "", STRATA)

slope_models <- readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/spline_models/",
                               p, "_full_model.rds"))[[sx]]

blups <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/spline_models/",
                           p, "_", sx, "_blups_full_model.txt"), 
                    sep = "\t", header = T, stringsAsFactors = F)

cluster_membership <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/",
                                        p, "_", sx, "_spline_full_models.txt"),
                                 sep = "\t", header = T, stringsAsFactors = F)

dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/indiv_qcd_data.rds")[[p]]
covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[[p]]

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220131_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)

age_at_diag_matrix <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_time_to_event_matrix.txt",
                                                   sep = "\t", header = T, stringsAsFactors = F)
colnames(age_at_diag_matrix) <- gsub("^X", "", colnames(age_at_diag_matrix))
age_at_diag_matrix <- age_at_diag_matrix[, 1:3]

TIME_VAR_COVARS <- c("data_provider")
TIME_INVAR_COVARS <- c("baseline_age", "age_sq",
                       "sex", "year_of_birth", "smoking_status")

# Wrangle data ----

# Add in age at first and last GP record to covariates
covars$eid <- as.character(covars$eid)
general_covars$eid <- as.character(general_covars$eid)
age_at_diag_matrix$eid <- as.character(age_at_diag_matrix$eid)

covars <- merge(covars, general_covars, by = "eid")
covars <- merge(covars, age_at_diag_matrix, by = "eid")

# Add in covariates to raw GP data
dat$eid <- as.character(dat$eid)
add_covs <- covars %>% select(any_of(c("eid", 
                                       TIME_VAR_COVARS, TIME_INVAR_COVARS)))
model_dat <- merge(dat, add_covs, by = "eid")
if (sx != "sex_comb") model_dat <- model_dat %>% filter(sex == sx)
model_dat <- model_dat %>% mutate(t = age_event - baseline_age)

# Count number of individuals in each cluster to make sampling quicker later
cluster_membership <- cluster_membership %>% 
  group_by(assigned_k) %>% 
  mutate(n_cluster = n())

# Plot raw data and model predictions for random individuals ----
## Get random IDs within each cluster ----

get_rand_eids <- function (n_each = 5) {
  # Return df of id and cluster number
  sampled_ids <- cluster_membership %>% 
    group_by(assigned_k) %>%
    sample_n(ifelse(n_cluster < n_each, n_cluster, n_each)) %>%
    select(eid, assigned_k)
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
    # Get covariates 
    select(any_of(c("eid", "max_t", 
                    TIME_VAR_COVARS, TIME_INVAR_COVARS))) %>% 
    distinct(across(all_of(c("eid", TIME_VAR_COVARS))), .keep_all = T) 
  
  # Timepoints to extend to
  ts <- lapply(1:nrow(new_data), FUN = function (i) { 
    seq(0, new_data$max_t[i], by = 0.25) })
  NT <- unlist(lapply(ts, function (x) length(x)))
  ts <- unlist(ts)
  
  new_data$ntimes <- NT
  new_data <- as.data.frame(lapply(new_data, rep, new_data$ntimes))
  new_data <- new_data %>% 
    mutate(t = ts,
           age_event = t + baseline_age,
           age_event_sq = age_event^2)
  
  # Predict new values
  fitted_results <- as.data.frame(predict(model_name,
                                          newdata = new_data))
  colnames(fitted_results) <- "fit"
  pred_df <- bind_cols(new_data, fitted_results)
  
  # At each time-point, average across time-varying covars
  pred_df <- pred_df %>% 
    group_by(eid, age_event, t) %>%
    summarise(fit = mean(fit)) 
  
  return (pred_df)
}

## Function to create plots ----

plot_predictions <- function (id_df) {
  raw_dat <- model_dat %>% filter(eid %in% id_df$eid)
  
  plot_dat <- create_prediction_df(slope_models, id_df$eid)
  plot_dat$cluster <- id_df$assigned_k[match(plot_dat$eid,
                                    id_df$eid)]
  plot_dat$cluster <- as.factor(as.character(plot_dat$cluster))
  
  # Order the ids by cluster number for plot
  id_levels <- id_df$eid[order(id_df$assigned_k)]
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
  highlight_ids <- ids_to_plot %>% group_by(assigned_k) %>% sample_n(5)
  highlight_ids <- unique(highlight_ids$eid)
  
  # Create fitted data
  plot_dat <- create_prediction_df(slope_models, ids_to_plot$eid)
  plot_dat$cluster <- ids_to_plot$assigned_k[match(plot_dat$eid,
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

# Plot cluster-level fitted trajectories ----

plot_fitted_traj_by_cluster <- function () {
  raw_dat_with_k <- model_dat %>%
    filter(eid %in% cluster_membership$eid) 
  raw_dat_with_k$cluster <- 
    cluster_membership$assigned_k[match(raw_dat_with_k$eid,
                                        cluster_membership$eid)]
  # Get nearest 0.25-yr time-bin
  summ_dat <- raw_dat_with_k %>% 
    mutate(t = plyr::round_any(t, 0.25, f = floor)) %>%
    filter(t <= 20) %>%
    group_by(cluster, t) %>% 
    summarise(mean_value = mean(value),
              sd_value = sd(value), 
              n = n()) %>%
    mutate(lci_mean = mean_value - 1.96*(sd_value/sqrt(n)),
           uci_mean = mean_value + 1.96*(sd_value/sqrt(n))) 
  
  # Get rolling average of these stats across 2 years 
  # because the plots are too choppy otherwise
  summ_dat <- summ_dat %>% 
    ungroup() %>% 
    group_by(cluster) %>% 
    arrange(t, .by_group = T) %>%
    mutate(interval_width = seq_along(t) - findInterval(t - 2, t),
           mean_value_rolled = rollapplyr(mean_value, interval_width, mean, 
                                          fill = NA),
           cluster = as.factor(as.character(cluster)))
  
  all_plot <- ggplot(summ_dat,
                     aes(x = t, y = mean_value, 
                         color = cluster, fill = cluster)) +
    # Faintly plot mean and SE in background
    geom_ribbon(aes(ymin = lci_mean, 
                    ymax = uci_mean), alpha = 0.2, 
                linetype = 0) +
    geom_line(aes(y = mean_value), alpha = 0.2) +
    # Add a thicker line for rolling average 
    geom_line(aes(y = mean_value_rolled), alpha = 1) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "Time from baseline measurement", 
         y = paste0("Mean and 95% C.I. of mean of ", p), 
         title = paste0("Observed trajectories in each cluster of phenotype: ",
                        p, " strata: ", sx))
  
  return (all_plot)
  
}

# Plot cluster-level observed trajectories ----

plot_mean_traj_by_cluster <- function () {
  raw_dat_with_k <- model_dat %>%
    filter(eid %in% cluster_membership$eid) 
  raw_dat_with_k$cluster <- 
    cluster_membership$assigned_k[match(raw_dat_with_k$eid,
                               cluster_membership$eid)]
  # Get nearest 0.25-yr age
  summ_dat <- raw_dat_with_k %>% 
    mutate(age_bin = plyr::round_any(age_event, 0.25, f = floor)) %>%
    group_by(cluster, age_bin) %>% 
    summarise(mean_value = mean(value),
              sd_value = sd(value), 
              n = n()) %>%
    mutate(lci_mean = mean_value - 1.96*(sd_value/sqrt(n)),
           uci_mean = mean_value + 1.96*(sd_value/sqrt(n))) 
  
  # Get rolling average of these stats across 2 years 
  # because the plots are too choppy otherwise
  summ_dat <- summ_dat %>% 
    ungroup() %>% 
    group_by(cluster) %>% 
    arrange(age_bin, .by_group = T) %>%
    mutate(interval_width = seq_along(age_bin) - 
             findInterval(age_bin - 2, age_bin),
           mean_value_rolled = rollapplyr(mean_value, interval_width, mean, 
                                        fill = NA),
           cluster = as.factor(as.character(cluster)))
  
  all_plot <- ggplot(summ_dat,
                     aes(x = age_bin, y = mean_value, 
                         color = cluster, fill = cluster)) +
    # Faintly plot mean and SE in background
    geom_ribbon(aes(ymin = lci_mean, 
                    ymax = uci_mean), alpha = 0.2, 
                linetype = 0) +
    geom_line(aes(y = mean_value), alpha = 0.2) +
    # Add a thicker line for rolling average 
    geom_line(aes(y = mean_value_rolled), alpha = 1) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "Age (years)", 
         y = paste0("Mean and 95% C.I. of mean of ", p), 
         title = paste0("Observed trajectories in each cluster of phenotype: ",
                        p, " strata: ", sx))
  
  return (all_plot)
  
}

# Plot cluster identity vs covariates ----

sex_col_palette <- c("#F8766D", "#00BFC4", "#C77CFF")
names(sex_col_palette) <- c("F", "M", "sex_comb")

plot_cluster_vs_covs <- function () {
  # Get covariates data
  covars_with_k <- covars %>%
    filter(eid %in% cluster_membership$eid) 
  covars_with_k$cluster <- 
    cluster_membership$assigned_k[match(covars_with_k$eid,
                               cluster_membership$eid)]
  covars_with_k$cluster <- as.factor(as.character(covars_with_k$cluster))
  
  # Number of men/women (individuals) in each cluster
  summ_sex <- covars_with_k %>% count(cluster, sex)
  sex_plot <- ggplot(summ_sex, 
                     aes(x = cluster, y = n, fill = sex)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(values = sex_col_palette) +
    labs(x = "Cluster", y = "Number of individuals") 
  
  # Year of birth stratified by cluster identity
  yob_plot <- ggplot(covars_with_k, 
                        aes(x = cluster, y = year_of_birth)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "Cluster", y = "Birth year") +
    theme(legend.position = "none")
  
  # Age at death stratified by cluster identity
  death_plot <- ggplot(covars_with_k, 
                     aes(x = cluster, y = age_at_death)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "Cluster", y = "Age at death (years)") +
    theme(legend.position = "none")
  
  # Age at first GP record stratified by cluster identity
  age_first_rec_plot <- ggplot(covars_with_k, 
                     aes(x = cluster, y = age_at_first_record)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "Cluster", y = "Age first seen in GP (years)") +
    theme(legend.position = "none")
  
  # Age at last GP record stratified by cluster identity
  age_last_rec_plot <- ggplot(covars_with_k, 
                     aes(x = cluster, y = age_at_last_record)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "Cluster", y = "Age last seen in GP (years)") +
    theme(legend.position = "none")
  
  # Baseline age stratified by cluster identity
  bl_age_plot <- ggplot(covars_with_k, 
                        aes(x = cluster, y = baseline_age)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "Cluster", y = paste0("Age ", p, 
                                            " first measured (years)")) +
    theme(legend.position = "none")
  
  # Baseline trait value
  bl_trait_plot <- ggplot(covars_with_k, 
                          aes(x = cluster, y = baseline_trait)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "Cluster", y = paste0("Baseline ", p)) +
    theme(legend.position = "none")
  
  # Number of follow-up measures
  max_y <- min(50, max(covars_with_k$FU_n))
  fu_n_plot <- ggplot(covars_with_k, 
                      aes(x = cluster, y = FU_n)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Dark2") + 
    scale_y_continuous(limits = c(0, max_y)) +
    labs(x = "Cluster", y = "# follow-up measures") +
    theme(legend.position = "none")
  
  # Number of follow-up years
  fu_yrs_plot <- ggplot(covars_with_k, 
                        aes(x = cluster, y = FUyrs)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "Cluster", y = "# years follow-up") +
    theme(legend.position = "none")
  
  # Smoking status in each cluster
  summ_smoking <- covars_with_k %>% 
    group_by(cluster, smoking_status) %>%
    summarise(n = n()) %>% 
    mutate(freq = n/sum(n))
  smoker_plot <- ggplot(summ_smoking, 
                     aes(x = cluster, y = freq, fill = smoking_status)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Cluster", y = "Proportion of individuals") 
  
  
  arranged_plots <- ggarrange(plotlist = list(sex_plot, yob_plot, death_plot,
                                              age_first_rec_plot, age_last_rec_plot,
                                              bl_age_plot, bl_trait_plot, 
                                              fu_n_plot, fu_yrs_plot,
                                              smoker_plot),
                              nrow = 2, ncol = 2)
  
  return (arranged_plots)
}

# Apply all plotting functions ----

# Modelled and observed trajectories of random ids
rand_plots <- plot_predictions(get_rand_eids(5))

# Modelled trajectories of random ids by cluster
cluster_sample_plots <- plot_samples_modelled()

# Refit trajectories within each cluster
refit_plots <- plot_fitted_traj_by_cluster()

# Modelled trajectories (all in cluster)
# modelled_all_plots <- plot_all_modelled(p, sx)

# Population-level observed trajectories
popn_plots <- plot_mean_traj_by_cluster()

# Cluster membership vs covariates 
cov_plots <- plot_cluster_vs_covs()

pdf(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/plots/trajectories/",
           p, "_", sx, "_trajectories.pdf"))
print(rand_plots)
print(cluster_sample_plots)
print(refit_plots)
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
k_info <- k_info[order(k_info$assigned_k), ]
row_annots <- data.frame(cluster = as.character(k_info$assigned_k))
rownames(row_annots) <- as.character(k_info$eid)

for_plot <- blups %>% filter(eid %in% k_info$eid)
rownames(for_plot) <- for_plot$eid
# Arrange by id and get rid of id column
for_plot <- for_plot[as.character(k_info$eid), -1]
for_plot <- as.matrix(for_plot)
# Replace parentheses text with blanks but retain (Intercept)
colnames(for_plot) <- gsub("\\s*\\([^\\)]+\\)", "", 
                           colnames(for_plot))

# For plotting, scale each column (RINT) to see differences better
for_plot <- apply(for_plot, 2, FUN = rint_x)

pdf(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/plots/assignment/",
           p, "_", sx, "_heatmaps.png"))
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

