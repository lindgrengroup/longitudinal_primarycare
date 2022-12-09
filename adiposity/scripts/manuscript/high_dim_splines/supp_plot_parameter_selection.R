# Author: Samvida S. Venkatesh
# Date: 24/05/2022

# Bits of code to plot parameter selection results

library(argparse)
library(splines)
library(cluster)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
theme_set(theme_bw())

RANDOM_SEED <- 160522
set.seed(RANDOM_SEED)

custom_teal_sequential <- c("#66BFBE", "#009593", "#005958")
names(custom_teal_sequential) <- c("2", "5", "10")
custom_four_diverge <- c("#D35C79", "#D9AB90", "#9FBCA4", "#009593")

main_filepath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/clustering/"

# Script to plot silhouette scores for combinations of K, L, M ----

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

# Plot silhouette plot given data frame of K, L, M, mean score, and S.D. score
plotSilScores <- function (df) {
  df <- df %>%
    mutate(K = as.numeric(K),
           L = as.factor(L),
           M = as.factor(M))
  resplot <- ggplot(df, aes(x = K, y = mean_silscore,
                        colour = L, fill = L)) +
    geom_point(aes(shape = M)) +
    geom_line(aes(group = interaction(L, M))) +
    geom_errorbar(aes(ymin = mean_silscore - 1.96*sd_silscore,
                      ymax = mean_silscore + 1.96*sd_silscore),
                  width = 0.1, alpha = 0.2) +
    scale_colour_manual(values = custom_teal_sequential) +
    scale_fill_manual(values = custom_teal_sequential) +
    scale_x_continuous(breaks = 2:8, labels = 2:8) +
    labs(y = "Mean (95% C.I.) across iterations silhouette score") +
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8))
  return (resplot)
}

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    all_combos <- list.files(paste0(main_filepath, p, "_", sx, "/parameter_selection/"),
                             pattern = "*.rds")
    centroid_dat <- lapply(all_combos, function (fname) {
      dat <- readRDS(paste0(main_filepath, p, "_", sx, "/parameter_selection/",
                                        fname))
      df_to_plot <- data.frame(K = dat$K, L = dat$L, M = dat$M,
                               mean_silscore = dat$silhouette_score$mean,
                               sd_silscore = dat$silhouette_score$sd)
      df_to_plot <- df_to_plot %>%
        mutate(K = as.numeric(K),
               L = as.character(L),
               M = as.character(M))
      return (df_to_plot)
    })
    all_parameters_dat <- bind_rows(centroid_dat)

    subset_dat <- all_parameters_dat %>% filter(K %in% c(3, 4, 5) &
                                                  L == "2")

    png(filename = paste0(main_filepath, 
                          "/silhouette_plots/zoomed_in_", p, "_", sx, ".png"),
        res = 300, units = "cm", height = 7, width = 7)
    print(plotSilScores(subset_dat))
    dev.off()

    # png(filename = paste0(main_filepath, "/silhouette_plots/", p, "_", sx, 
    #                       "_silhouette_plot.png"),
    #     res = 300, units = "cm", height = 7, width = 7)
    # print(plotSilScores(all_parameters_dat))
    # dev.off()
  })
})

# Script to plot combined iteration centroids for all values of L, M, when given K ----

## Read data ----

parser <- ArgumentParser()
parser$add_argument("--phenotype", required = TRUE,
                    help = "Phenotype to model")
parser$add_argument("--ss", required = TRUE,
                    help = "Sex strata")
args <- parser$parse_args()

PHENO <- args$phenotype
SEX_STRATA <- args$ss
K <- 2:8

centroid_dat <- lapply(K, function (k) {
  all_combos <- list.files(paste0(main_filepath, PHENO, "_", SEX_STRATA, 
                                  "/parameter_selection/"),
                           pattern = paste0("^K", k))
  res_list <- lapply(all_combos, function (fname) {
    dat <- readRDS(paste0(main_filepath, PHENO, "_", SEX_STRATA, 
                          "/parameter_selection/",
                          fname))
    metadat <- data.frame(L = as.character(dat$L), 
                          M = as.character(dat$M),
                          K = as.character(k))
    dat <- dat$cluster_centroids
    return (list(dat = dat,
                 metadat = metadat))
  })
  return (res_list)
})
names(centroid_dat) <- paste0("k", K)

# Modelling results for distance calculations 
model_dat <- readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/results/with_rvar_fit_objects_", 
                            PHENO, "_", SEX_STRATA, ".rds"))

## Calculate mean matrix of coefficients needed for distance calculations ----

B <- model_dat$B
spline_posteriors <- model_dat$spline_posteriors
model_resid_var <- model_dat$resid_var

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

per_k_plots <- lapply(K, function (k) {
  sanity_check_cluster_centres <- lapply(centroid_dat[[paste0("k", k)]], function (dlist) {
    dat_plot <- getPredValuesClusterCentroid(dlist$dat[, -1])
    dat_plot$parameter <- paste0("K", k, 
                                 "_L", dlist$metadat$L, 
                                 "_M", dlist$metadat$M)
    return (dat_plot)
  })
  sanity_check_cluster_centres <- bind_rows(sanity_check_cluster_centres) %>%
    mutate(parameter = factor(parameter, 
                              levels = paste0("K", k, "_",
                                              c("L2_Mrandom", "L2_M1", "L2_M2", "L2_M5", "L2_M10",
                                                "L5_Mrandom", "L5_M1", "L5_M2", "L5_M5", "L5_M10",
                                                "L10_Mrandom", "L10_M1", "L10_M2", "L10_M5", "L10_M10"))))
  
  usecolpal <- colorRampPalette(custom_four_diverge)(k)
  names(usecolpal) <- 1:k
  
  clust_centres_plot <- ggplot(sanity_check_cluster_centres, 
                               aes(x = t_diff, y = pred_value, 
                                   col = clust, fill = clust)) +
    facet_wrap(~parameter, ncol = 5) +
    geom_line() +
    scale_color_manual(values = usecolpal, guide = "none") +
    scale_fill_manual(values = usecolpal, guide = "none") +
    labs(x = "Days from first measurement", 
         y = "Cluster centroid predicted value") + 
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8))
  
  png(filename = paste0(main_filepath, "cluster_centroid_plots/", 
                        PHENO, "_", SEX_STRATA, "_K", k, "_cluster_centroids.png"),
      res = 300, units = "cm", height = 20, width = 15)
  print(clust_centres_plot)
  dev.off()
  
})

# Script to plot combined iteration centroids for matched values of K,L,M old vs new initialisation methods ----

## Read data ----

parser <- ArgumentParser()
parser$add_argument("--phenotype", required = TRUE,
                    help = "Phenotype to model")
parser$add_argument("--ss", required = TRUE,
                    help = "Sex strata")
args <- parser$parse_args()

PHENO <- args$phenotype
SEX_STRATA <- args$ss
K <- 2:8

centroid_dat <- lapply(K, function (k) {
  og_combos <- list.files(paste0(main_filepath, PHENO, "_", SEX_STRATA, 
                                  "/parameter_selection/"),
                           pattern = paste0("^K", k))
  og_list <- lapply(og_combos, function (fname) {
    dat <- readRDS(paste0(main_filepath, PHENO, "_", SEX_STRATA, 
                          "/parameter_selection/",
                          fname))
    metadat <- data.frame(L = as.character(dat$L), 
                          M = as.character(dat$M),
                          K = as.character(k),
                          type = "old")
    dat <- dat$cluster_centroids
    return (list(dat = dat,
                 metadat = metadat))
  })
  
  new_combos <- list.files(paste0(main_filepath, PHENO, "_", SEX_STRATA, 
                                 "/parameter_selection/medoid_initialisation/"),
                          pattern = paste0("^K", k))
  new_list <- lapply(new_combos, function (fname) {
    dat <- readRDS(paste0(main_filepath, PHENO, "_", SEX_STRATA, 
                          "/parameter_selection/medoid_initialisation/",
                          fname))
    metadat <- data.frame(L = as.character(dat$L), 
                          M = as.character(dat$M),
                          K = as.character(k),
                          type = "new")
    dat <- dat$cluster_centroids
    return (list(dat = dat,
                 metadat = metadat))
  })
  return (list(old_centroids = og_list,
               new_centroids = new_list))
})
names(centroid_dat) <- paste0("k", K)

# Modelling results for distance calculations 
model_dat <- readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/results/fit_objects_", 
                            PHENO, "_", SEX_STRATA, ".rds"))

## Calculate mean matrix of coefficients needed for distance calculations ----

B <- model_dat$B
spline_posteriors <- model_dat$spline_posteriors
model_resid_var <- model_dat$resid_var

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

per_k_plots <- lapply(K, function (k) {
  old_cluster_centres <- lapply(centroid_dat[[paste0("k", k)]][["old_centroids"]], function (dlist) {
    dat_plot <- getPredValuesClusterCentroid(dlist$dat[, -1])
    dat_plot$parameter <- paste0("K", k, 
                                 "_L", dlist$metadat$L, 
                                 "_M", dlist$metadat$M)
    dat_plot$type <- dlist$metadat$type
    return (dat_plot)
  })
  old_cluster_centres <- bind_rows(old_cluster_centres)
  
  new_cluster_centres <- lapply(centroid_dat[[paste0("k", k)]][["new_centroids"]], function (dlist) {
    dat_plot <- getPredValuesClusterCentroid(dlist$dat[, -1])
    dat_plot$parameter <- paste0("K", k, 
                                 "_L", dlist$metadat$L, 
                                 "_M", dlist$metadat$M)
    dat_plot$type <- dlist$metadat$type
    return (dat_plot)
  })
  new_cluster_centres <- bind_rows(new_cluster_centres)
  
  for_plot <- bind_rows(old_cluster_centres, new_cluster_centres) %>%
    mutate(parameter = factor(parameter, 
                              levels = paste0("K", k, "_",
                                              c("L2_Mrandom", "L2_M1", "L2_M2", "L2_M5", "L2_M10",
                                                "L5_Mrandom", "L5_M1", "L5_M2", "L5_M5", "L5_M10",
                                                "L10_Mrandom", "L10_M1", "L10_M2", "L10_M5", "L10_M10"))),
           type = factor(type))
  
  usecolpal <- colorRampPalette(custom_four_diverge)(k)
  names(usecolpal) <- 1:k
  
  clust_centres_plot <- ggplot(for_plot, 
                               aes(x = t_diff, y = pred_value, 
                                   col = clust, fill = clust,
                                   linetype = type)) +
    facet_wrap(~parameter, ncol = 5) +
    geom_line() +
    scale_color_manual(values = usecolpal, guide = "none") +
    scale_fill_manual(values = usecolpal, guide = "none") +
    scale_linetype_manual(values = c("new" = 1, "old" = 2)) +
    labs(x = "Days from first measurement", 
         y = "Cluster centroid predicted value") + 
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8))
  
  png(filename = paste0(main_filepath, "cluster_centroid_plots/", 
                        PHENO, "_", SEX_STRATA, "_K", k, "_cluster_centroids.png"),
      res = 300, units = "cm", height = 20, width = 15)
  print(clust_centres_plot)
  dev.off()
  
})


