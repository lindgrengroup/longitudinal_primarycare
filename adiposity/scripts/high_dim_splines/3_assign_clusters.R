# Author: Samvida S. Venkatesh
# Date: 24/05/2022

# Script to plot silhouette scores for combinations of K, L, M ----
# 
# library(tidyverse)
# theme_set(theme_bw())
# 
# RANDOM_SEED <- 160522
# set.seed(RANDOM_SEED)
# 
# custom_teal_sequential <- c("#66BFBE", "#009593", "#005958")
# names(custom_teal_sequential) <- c("2", "5", "10")
# 
# PHENOTYPES <- c("BMI", "Weight")
# SEX_STRATA <- c("F", "M", "sex_comb")
# 
# main_filepath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/"
# 
# # Plot silhouette plot given data frame of K, L, M, mean score, and S.D. score
# plotSilScores <- function (df) {
#   df <- df %>% 
#     mutate(K = as.numeric(K),
#            L = as.factor(L),
#            M = as.factor(M))
#   resplot <- ggplot(df, aes(x = K, y = mean_silscore,
#                         colour = L, fill = L)) +
#     geom_point(aes(shape = M)) +
#     geom_line(aes(group = interaction(L, M))) +
#     geom_errorbar(aes(ymin = mean_silscore - 1.96*sd_silscore,
#                       ymax = mean_silscore + 1.96*sd_silscore),
#                   width = 0.1, alpha = 0.2) +
#     scale_colour_manual(values = custom_teal_sequential) +
#     scale_fill_manual(values = custom_teal_sequential) +
#     scale_x_continuous(breaks = 2:8, labels = 2:8) +
#     labs(y = "Mean (95% C.I.) across iterations silhouette score") +
#     theme(axis.text = element_text(size = 6),
#           axis.title = element_text(size = 8),
#           legend.text = element_text(size = 6),
#           legend.title = element_text(size = 8))
#   return (resplot)
# }
# 
# lapply(PHENOTYPES, function (p) {
#   lapply(SEX_STRATA, function (sx) {
#     all_combos <- list.files(paste0(main_filepath, p, "_", sx, "/parameter_selection/"),
#                              pattern = "*.rds")
#     centroid_dat <- lapply(all_combos, function (fname) {
#       dat <- readRDS(paste0(main_filepath, p, "_", sx, "/parameter_selection/",
#                                         fname))
#       df_to_plot <- data.frame(K = dat$K, L = dat$L, M = dat$M,
#                                mean_silscore = dat$silhouette_score$mean,
#                                sd_silscore = dat$silhouette_score$sd)
#       df_to_plot <- df_to_plot %>%
#         mutate(K = as.numeric(K),
#                L = as.character(L),
#                M = as.character(M))
#       return (df_to_plot)
#     })
#     all_parameters_dat <- bind_rows(centroid_dat)
#     
#     # subset_dat <- all_parameters_dat %>% filter(K %in% c(3, 4, 5) &
#     #                                               L == "2")
#     # 
#     # png(filename = paste0(main_filepath, p, "_", sx, "/parameter_selection/zoomed_in_silhouette_plot.png"),
#     #     res = 300, units = "cm", height = 7, width = 7)
#     # print(plotSilScores(subset_dat))
#     # dev.off()
#     
#     png(filename = paste0(main_filepath, p, "_", sx, "/parameter_selection/silhouette_plot.png"),
#         res = 300, units = "cm", height = 7, width = 7)
#     print(plotSilScores(all_parameters_dat))
#     dev.off()
#   })
# })

# Main script ----

library(argparse)
library(splines)
library(zoo)
library(tidyverse)
theme_set(theme_bw())

RANDOM_SEED <- 160522
set.seed(RANDOM_SEED)

custom_four_diverge <- c("#D35C79", "#D9AB90", "#9FBCA4", "#009593")
names(custom_four_diverge) <- c("1", "2", "3", "4")

parser <- ArgumentParser()
parser$add_argument("--phenotype", required = TRUE,
                    help = "Phenotype to model")
parser$add_argument("--sex_strata", required = TRUE,
                    help = "Sex strata")
parser$add_argument("--K", required = TRUE,
                    help = "Chosen number of clusters")
parser$add_argument("--L", required = TRUE,
                    help = "Chosen minimum L")
parser$add_argument("--M", required = TRUE,
                    help = "Chosen initialisation strategy")
args <- parser$parse_args()

PHENO <- args$phenotype
SEX_STRATA <- args$sex_strata
K <- as.numeric(args$K)
L <- args$L
M <- args$M

main_filepath <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/",
                        PHENO, "_", SEX_STRATA, "/")

# Read data ----

# Centroids from chosen clustering method
clust_centroids <- readRDS(paste0(main_filepath, 
                                  "parameter_selection/K", K, "_L", L, "_M", M, ".rds"))
clust_centroids <- clust_centroids$cluster_centroids

# Modelling results for distance calculations 
model_dat <- readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/results/fit_objects_", 
                            PHENO, "_", SEX_STRATA, ".rds"))

# Covariates to compare between training and validation sets
covars_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[[PHENO]] %>%
  mutate(eid = as.character(eid))
general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F) %>%
  mutate(eid = as.character(eid))

# Assignment of ids to training or validation sets 
training_ids <- read.table(paste0(main_filepath, "/ids_training.txt"),
                           sep = "\t", header = F, stringsAsFactors = F)$V1
training_ids <- data.frame(eid = as.character(training_ids),
                           id_type = "training")
validation_ids <- read.table(paste0(main_filepath, "/ids_validation.txt"),
                           sep = "\t", header = F, stringsAsFactors = F)$V1
validation_ids <- data.frame(eid = as.character(validation_ids),
                           id_type = "validation")
id_classification <- bind_rows(training_ids, validation_ids)

# Original population-level data
orig_popn_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/indiv_qcd_data.rds")[[PHENO]] %>%
  mutate(eid = as.character(eid))

# Calculate mean matrix of coefficients needed for distance calculations ----

B <- model_dat$B
spline_posteriors <- model_dat$spline_posteriors
model_resid_var <- model_dat$resid_var

# Create mean coefficient matrix
# Matrix of means (id x basis)
mn_mat <- lapply(spline_posteriors, function (spobj) {
  return (as.data.frame(t(spobj$mu)))
})
mn_mat <- bind_rows(mn_mat)
rownames(mn_mat) <- names(spline_posteriors)

# Baseline data before clustering (subtract intercept - value at t0)
mn_mat <- mn_mat - mn_mat[, 1]

# Plot final cluster centroids to sanity check ----

## Plotting functions ----

## Predictions
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
    mutate(t_diff = as.numeric(t_diff))
  return (for_plot)
}

## Apply ----

sanity_check_cluster_centres <- getPredValuesClusterCentroid(clust_centroids[, -1])

clust_centres_plot <- ggplot(sanity_check_cluster_centres, 
                             aes(x = t_diff, y = pred_value, 
                                 col = clust, fill = clust)) +
  geom_line() +
  scale_color_manual(values = custom_four_diverge, guide = "none") +
  scale_fill_manual(values = custom_four_diverge, guide = "none") +
  labs(x = "Days from first measurement", 
       y = "Cluster centroid predicted value") + 
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8))
  
png(filename = paste0(main_filepath, "/plots/cluster_centroids.png"),
    res = 300, units = "cm", height = 7, width = 7)
print(clust_centres_plot)
dev.off()

# Functions to assign all individuals to their closest cluster ----

# Given a vector and a matrix, get row index in matrix closest to the vector
getClosestMedoid <- function (x1, mat_xs) {
  x1 <- unlist(x1)
  all_dists <- t(apply(mat_xs, 1, function (x) x - x1))^2
  return (which.min(rowSums(all_dists)))
}

# Given a list of ids and centroid positions, assign ids to closest medoid
getClusterBelonging <- function (id_list, centroid_mat) {
  id_means <- mn_mat[id_list, ]
  
  # Euclidean distance from each id to centroid
  closest_medoid <- lapply(1:nrow(id_means), 
                           function (ii) getClosestMedoid(id_means[ii, ], centroid_mat))
  
  res <- data.frame(eid = rownames(id_means),
                    clust = unlist(closest_medoid))
  
  return (res)
}

final_clust_assignments <- getClusterBelonging(id_list = id_classification$eid,
                                               centroid_mat = clust_centroids[, -1])

# Compare between training and validation sets ----

# 1. Fraction of individuals assigned to each cluster
# 2. Covariate distributions
# 3. Mean population level trajectories

full_dat <- left_join(final_clust_assignments, id_classification, by = "eid")
full_dat <- left_join(full_dat, covars_dat, by = "eid")
full_dat <- left_join(full_dat, general_covars, by = "eid")

sumstats_table <- full_dat %>% 
  group_by(clust, id_type) %>%
  summarise(count = n(), 
            mean_FU_n = signif(mean(FU_n), 3), 
            sd_FU_n = signif(sd(FU_n), 3),
            mean_FUyrs = signif(mean(FUyrs), 3), 
            sd_FUyrs = signif(sd(FUyrs), 3),
            mean_bl_age = signif(mean(baseline_age), 3), 
            sd_bl_age = signif(sd(baseline_age), 3),
            median_bl_trait = signif(median(baseline_trait), 3), 
            iqr_bl_trait = paste(signif(quantile(baseline_trait, 0.25), 3), 
                                 signif(quantile(baseline_trait, 0.75), 3),
                                 sep = ", ")) 

write.table(sumstats_table,
            paste0(main_filepath, "/training_vs_validation_cluster_properties.txt"),
            sep = "\t", row.names = F, quote = F)

## Plotting functions ----

to_plot <- full_dat %>%
  mutate(clust = as.factor(as.character(clust)),
         id_type = as.factor(id_type))

violinCovarPlot <- function (dat, covar_name) {
  resplot <- ggplot(dat, aes(x = clust, y = !!as.symbol(covar_name),
                             group = interaction(clust, id_type))) +
    geom_violin(aes(fill = clust, 
                    alpha = id_type), 
                position = position_dodge(1)) +
    geom_boxplot(position = position_dodge(1), width = 0.1) + 
    scale_fill_manual(values = custom_four_diverge, guide = "none") +
    scale_alpha_manual(values = c(0.5, 1), guide = "none") + 
    labs(x = "Cluster", y = covar_name) + 
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8))
  
  return (resplot)
}

# Apply covariate plotting to all covariates 
lapply(c("baseline_age", "baseline_trait", "FU_n", "FUyrs"), function (cvname) {
  
  png(filename = paste0(main_filepath, "/plots/training_vs_validation_",
                        cvname, ".png"),
      res = 300, units = "cm", height = 7, width = 7)
  print(violinCovarPlot(to_plot, cvname))
  dev.off()

})

to_plot <- left_join(orig_popn_dat, final_clust_assignments, by = "eid")
to_plot <- left_join(to_plot, id_classification, by = "eid")
to_plot <- to_plot %>%
  filter(!is.na(clust) & !is.na(id_type)) %>%
  mutate(clust = as.factor(as.character(clust)),
         id_lty = ifelse(id_type == "validation", 1, 2)) 

plotPopnTrajCluster <- function (dat) {
  
  # Round ages to nearest 0.25 yrs
  summ_dat <- dat %>% 
    mutate(age_bin = plyr::round_any(age_event, 0.25, f = round)) %>%
    filter(age_bin >= 30 & age_bin <= 70) %>%
    group_by(clust, age_bin, id_lty) %>% 
    summarise(mean_value = mean(value))
  
  # Get rolling average mean across 2 years 
  summ_dat <- summ_dat %>% 
    ungroup() %>% 
    group_by(clust, id_lty) %>% 
    arrange(age_bin, .by_group = T) %>%
    mutate(interval_width = seq_along(age_bin) - 
             findInterval(age_bin - 2, age_bin),
           mean_value_rolled = rollapply(mean_value, interval_width, mean, 
                                         fill = NA),
           clust = factor(as.character(clust)),
           id_lty = factor(as.character(id_lty)))
  
  all_plot <- ggplot(summ_dat,
                     aes(x = age_bin, y = mean_value_rolled, 
                         color = clust, linetype = id_lty)) +
    # Add a thick line for rolling average 
    geom_line() +
    scale_color_manual(values = custom_four_diverge) +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE))
  
  return (all_plot)
}

png(filename = paste0(main_filepath, "/plots/training_vs_validation_trajectories.png"),
    res = 300, units = "cm", height = 7, width = 7)
print(plotPopnTrajCluster(to_plot))
dev.off()

# Write file for GWAS ----

# IDs that passed sample QC
ids_passed_qc <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/sample_qc/", 
                             PHENO, "_", SEX_STRATA, "_ids_passed_qc.txt"),
                      sep = "\t", header = T)
ids_passed_qc$IID <- as.character(ids_passed_qc$IID)

PCs <- paste0("PC", c(1:21))
COVARS_LIST <- c("UKB_assmt_centre", "genotyping.array", 
                 "year_of_birth", "baseline_age", "age_sq", 
                 "FU_n", "FUyrs", "sex", PCs)

# Wrangle data to wide form (k1, k2, etc. for cluster belonging)
to_write <- full_dat %>% select(all_of(c("eid", "clust", COVARS_LIST))) %>%
  filter(eid %in% ids_passed_qc$IID) %>%
  mutate(clust = paste0("k", clust)) %>%
  pivot_wider(id_cols = -clust, 
              names_from = clust, 
              values_from = clust, 
              values_fn = function (x) 1, 
              values_fill = 0)

# Create group columns
to_write <- to_write %>% 
  mutate(k1_k2 = ifelse(k1 + k2 > 0, 1, 0))

# Write results to table
write.table(to_write,
            paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/cluster_membership_", 
                   PHENO, "_", SEX_STRATA, ".txt"),
            sep = " ", row.names = F, quote = F)
