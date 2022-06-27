# Author: Samvida S. Venkatesh
# Adapted from: George Nicholson
# Date: 16/05/22

library(argparse)
library(splines)
library(cluster)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
theme_set(theme_bw())

RANDOM_SEED <- 160522
set.seed(RANDOM_SEED)

parser <- ArgumentParser()
parser$add_argument("--phenotype", required = TRUE,
                    help = "Phenotype to model")
parser$add_argument("--ss", required = TRUE,
                    help = "Sex strata")
parser$add_argument("--nclust", required = TRUE,
                    default = 5,
                    help = "Number of clusters")
parser$add_argument("--L", required = TRUE,
                    help = "Only sample from individuals with at least L measurements")
parser$add_argument("--M", required = TRUE,
                    help = "Initialise clusters with diff. at M years post-baseline")

args <- parser$parse_args()

PHENO <- args$phenotype
SEX_STRATA <- args$ss
NCLUST <- as.numeric(args$nclust)

# For sampling scheme
NSAMPLES <- 10000 # number of individuals to sample each iteration
S <- 9 # number of times to draw samples from population
L <- as.numeric(args$L) # only sample from individuals with at least L measurements

# For clustering
M <- args$M
if (M != "random") {
  M <- as.numeric(args$M) # initialise clusters with n-tile diff. at M years post-baseline
}

plotdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/plots/sample_clustering_scheme"
resdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/"

# Load data ----

model_dat <- readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/results/fit_objects_", 
                            PHENO, "_", SEX_STRATA, ".rds"))
covars_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[[PHENO]]

# Get data relevant for sampling scheme ----

B <- model_dat$B
spline_posteriors <- model_dat$spline_posteriors
model_resid_var <- model_dat$resid_var

# Create mean and S.D. matrix
# Matrix of means (id x basis)
mn_mat <- lapply(spline_posteriors, function (spobj) {
  return (as.data.frame(t(spobj$mu)))
})
mn_mat <- bind_rows(mn_mat)
rownames(mn_mat) <- names(spline_posteriors)

# Matrix of standard deviations (id x basis)
sd_mat <- lapply(spline_posteriors, function (spobj) {
  return (sqrt(diag(spobj$Sig * model_resid_var)))
})
sd_mat <- bind_rows(sd_mat)
sd_mat <- as.data.frame(sd_mat)
rownames(sd_mat) <- names(spline_posteriors)

# Minimum number of measurements
valid_ids <- as.character(covars_dat$eid[covars_dat$FU_n >= L])
valid_ids <- rownames(model_dat$mn_mat)[rownames(model_dat$mn_mat) %in% valid_ids]

mn_mat <- model_dat$mn_mat[valid_ids, ]
sd_mat <- model_dat$sd_mat[valid_ids, ]
spline_posteriors <- spline_posteriors[valid_ids]

# Baseline data before clustering (subtract intercept - value at t0)
mn_mat <- mn_mat - mn_mat[, 1]

# Function for distance matrix calculation ----

# Given vector of ids, calculate distances between each pair of ids in vector
getCustDistMat <- function (id_list) {
  cat("######## Building distance matrix", "\n")
  sub_mn <- mn_mat[id_list, ]
  sub_sd <- sd_mat[id_list, ]
  # Initialise empty distance matrix
  dist_res <- matrix(NA, length(id_list), length(id_list),
                     dimnames = list(id_list, id_list))
  # Loop through each id
  for (j in 1:length(id_list)) {
    mn_curr <- unlist(sub_mn[j, ])
    mn_diff <- t(apply(sub_mn[j:length(id_list), ], 1, function (x) x - mn_curr))^2
    
    v_curr <- unlist(sub_sd[j, ]^2)
    v_all <- t(apply(sub_sd[j:length(id_list), ], 1, function (x) x^2 + v_curr))
    
    dist_curr <- sqrt(rowSums(mn_diff / v_all))
    dist_res[j, j:length(id_list)] <- dist_res[j:length(id_list), j] <- dist_curr
  }
  # Return matrix
  return (dist_res)
}

# Functions for clustering ----

# Given a list of ids, calculate which id is closest to the group centre
getMedoid <- function (id_list) {
  id_means <- mn_mat[id_list, ]
  gp_centroid <- colMeans(id_means)
  # Euclidean distance from each id to centroid
  mn_diff <- t(apply(id_means, 1, function (x) x - gp_centroid))^2
  medoid_id <- id_list[which.min(rowSums(mn_diff))]
  return (medoid_id)
}

# Given a vector and a matrix, get row index in matrix closest to the vector
getClosestMedoid <- function (x1, mat_xs) {
  x1 <- unlist(x1)
  all_dists <- t(apply(mat_xs, 1, function (x) x - x1))^2
  return (which.min(rowSums(all_dists)))
}

# Given a list of ids and medoids, assign ids to closest medoid
getClusterBelonging <- function (id_list, medoid_list) {
  id_means <- mn_mat[id_list, ]
  medoid_pos <- mn_mat[medoid_list, ]
  
  # Euclidean distance from each id to centroid
  closest_medoid <- lapply(1:nrow(id_means), 
                          function (ii) getClosestMedoid(id_means[ii, ], medoid_pos))
  
  res <- data.frame(eid = rownames(id_means),
                    init_clust = unlist(closest_medoid))
  
  return (res)
}

## Calculate cluster medoids for initialisation, given list of individuals
getInitialCentres <- function (id_list, myrs = M,
                               ncentres = NCLUST) {
  
  cat("######## Initialising cluster centres", "\n")
  
  if (M == "random") {
    
    # Sample initial medoids randomly from list
    init_meds <- sample(id_list, NCLUST, replace = F)
    index_ret <- match(init_meds, id_list)
    
    # Assign remaining ids to closest initial centre
    og_assignment <- getClusterBelonging(id_list, init_meds)
    
  } else {
    # Predict full trajectory for all individuals (returns id x time matrix)
    pred_mat <- t(apply(mn_mat[id_list, ], 1, 
                        function (x) model_dat$B %*% x))
    
    # Get fold-change between 1st and M*365th measurement for each individual
    # and assign to an n-tile
    if (ncol(pred_mat) < myrs*365) stop("Set smaller M")
    diffs_assess <- (pred_mat[, myrs*365] - pred_mat[, 1]) / pred_mat[, 1]
    get_id_q <- data.frame(eid = rownames(pred_mat),
                           fc = diffs_assess) %>%
      mutate(init_clust = ntile(fc, ncentres))
    
    init_meds <- get_id_q %>% group_by(init_clust) %>%
      summarise(gp_medoid = getMedoid(eid))
    
    # Return indices of medoids and initial cluster assignment
    og_assignment <- get_id_q[, c("eid", "init_clust")]
    index_ret <- match(init_meds$gp_medoid, id_list)
    
  }
  
  return (list(og_assignment = og_assignment,
               medoid_index = index_ret))
}

## Perform PAM clustering (partitioning around medoids) and return clusters

performClustering <- function (dist_mat, ncentres = NCLUST, init_centres) {
  clust_res <- pam(dist_mat, k = ncentres, diss = T, 
                   medoids = init_centres, 
                   keep.diss = F, cluster.only = T)
  
  res <- data.frame(eid = names(clust_res),
                    final_clust = clust_res)
  cat("######## Clustering done!", "\n", "\n")
  return (res)
}

# Plotting functions ----

## Cluster centroid trajectories ----

## Predictions
getPredValuesClusterCentroid <- function (coef_mat) {
  pred_vals <- t(apply(coef_mat, 1, 
                       function (x) model_dat$B %*% x))
  # Wrangle into ggplot format
  for_plot <- as.data.frame(pred_vals)
  colnames(for_plot) <- paste0("d", 1:ncol(for_plot))
  for_plot$final_clust <- factor(as.character(1:NCLUST))
  
  for_plot <- for_plot %>% pivot_longer(cols = -final_clust,
                                        names_to = "t_diff", 
                                        names_prefix = "d", 
                                        values_to = "pred_value") %>%
    mutate(t_diff = as.numeric(t_diff))
  return (for_plot)
}

## Cluster centroid mean and S.D., wrangled into long format for plot
getPlotDatClusterCentroids <- function (grouped_id_list) {
  
  full_dat <- mn_mat[grouped_id_list$eid, ]
  full_dat$eid <- rownames(full_dat)
  full_dat$final_clust <- 
    grouped_id_list$final_clust[match(full_dat$eid, grouped_id_list$eid)]
  
  centroid_means <- full_dat %>% group_by(final_clust) %>%
    summarise(across(-eid, mean))
  centroid_sds <- full_dat %>% group_by(final_clust) %>%
    summarise(across(-eid, sd))
  
  # Predict full trajectory for centroids (returns centroid x time matrix)
  pred_mns <- getPredValuesClusterCentroid(centroid_means[, -1]) %>%
    rename(pred_mean_value = pred_value)
  pred_locis <- getPredValuesClusterCentroid(centroid_means[, -1] - 1.96*centroid_sds[, -1]) %>%
    rename(pred_loci_value = pred_value)
  pred_upcis <- getPredValuesClusterCentroid(centroid_means[, -1] + 1.96*centroid_sds[, -1]) %>%
    rename(pred_upci_value = pred_value)
  
  # Wrangle into ggplot format
  for_plot <- full_join(pred_mns, pred_locis, by = c("final_clust", "t_diff"))
  for_plot <- full_join(for_plot, pred_upcis, by = c("final_clust", "t_diff"))
  
  return (for_plot)
}

## Confusion matrix plots for initial vs final cluster assignment ----

getPlotDatConfusionMat <- function (cmat) {
  cmat <- cmat %>% mutate(init_clust = factor(init_clust, 
                                                 levels = 1:NCLUST),
                          final_clust = factor(final_clust, 
                                                  levels = NCLUST:1))
  res <- ggplot(cmat, aes(x = init_clust, y = final_clust, fill = count)) +
    geom_tile() + 
    geom_text(aes(label = count), size = 2) +
    scale_fill_gradient(low = "white", high = "#009194", guide = "none") +
    labs(x = "Initial", y = "Final") 
  return (res)
}

# Apply clustering S times ----

itered_clustering <- lapply(1:S, function (si) {
  cat(paste0("Running iteration #", si), "\n")
  # Get ids
  id_list <- sample(valid_ids, NSAMPLES, replace = F)
  
  # Initialisation
  init_centres <- getInitialCentres(id_list, myrs = M, ncentres = NCLUST)
  # Cluster
  clust_res <- performClustering(dist_mat = getCustDistMat(id_list), 
                                 ncentres = NCLUST, 
                                 init_centres = init_centres$medoid_index)
  
  # Add in initial clustering assignment
  res <- left_join(clust_res, init_centres$og_assignment, by = "eid")
  return (res)
})
names(itered_clustering) <- paste0("iter", 1:S)

saveRDS(itered_clustering, 
        paste0(resdir, "sample_clustering_scheme_", PHENO, "_", SEX_STRATA, 
               "_L", L, "_M", M, 
               ".rds"))

# Plot results ----

## Centroid mean trajectories ---
centroid_traj <- lapply(1:S, function (si) {
  plot_res <- getPlotDatClusterCentroids(itered_clustering[[paste0("iter", si)]])
  plot_res$iteration <- as.character(si)
  return (plot_res)
})
centroid_traj <- bind_rows(centroid_traj)

centroid_traj_plot <- ggplot(centroid_traj, 
                             aes(x = t_diff, y = pred_mean_value,
                                 col = final_clust)) +
  facet_wrap(~iteration) +
  geom_line() +
  geom_ribbon(aes(ymin = pred_loci_value, ymax = pred_upci_value,
                  fill = final_clust),
              alpha = 0.1, linetype = 0) +
  scale_color_brewer(palette = "Set1", guide = "none") +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  labs(x = "Days from first measurement", 
       y = "Cluster centroid predicted value")

## Confusion matrices for initial vs final clustering assignment ----

confusion_mat_plots <- lapply(1:S, function (si) {
  cmat <- as.data.frame(table(itered_clustering[[si]]$init_clust, 
                              itered_clustering[[si]]$final_clust))
  colnames(cmat) <- c("init_clust", "final_clust", "count")
  resplot <- getPlotDatConfusionMat(cmat)
  return (resplot)
})

# Print all plots ----

pdf(paste0(plotdir, "sample_clustering_scheme_", PHENO, "_", SEX_STRATA, "_L", L, "_M", M, 
           ".pdf"), onefile = T)
print(centroid_traj_plot)
# Arrange confusion matrix plots
print(ggarrange(plotlist = confusion_mat_plots, ncol = 3, nrow = 3))
dev.off()
