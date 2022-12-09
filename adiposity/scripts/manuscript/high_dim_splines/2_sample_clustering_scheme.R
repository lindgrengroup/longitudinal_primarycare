# Author: Samvida S. Venkatesh
# Adapted from: George Nicholson
# Date: 16/05/22

# ONLY RUN THE FOLLOWING ONCE TO SPLIT DATA INTO 80% CLUSTER DISCOVERY SET
# AND 20% CLUSTER VALIDATION SET

# Split data into training and validation ----

# library(tidyverse)
# RANDOM_SEED <- 160522
# set.seed(RANDOM_SEED)
# 
# PHENOTYPES <- c("BMI", "Weight")
# SEX_STRATA <- c("F", "M", "sex_comb")
# 
# main_filepath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/"
# raw_dat <- readRDS(paste0(main_filepath, "data/dat_to_model.rds"))
# 
# lapply(PHENOTYPES, function (p) {
#   lapply(SEX_STRATA, function (sx) {
#     df <- raw_dat[[p]][[sx]] 
#     
#     ids_to_split <- unique(df$eid)
#     NTRAIN <- round(0.8*length(ids_to_split))
#     
#     training_ids <- sample(ids_to_split, NTRAIN, replace = F)
#     validation_ids <- ids_to_split[!ids_to_split %in% training_ids]
#     
#     write.table(training_ids, 
#                 paste0(main_filepath, "clustering/",
#                        p, "_", sx, "/ids_training.txt"),
#                 sep = "\t", row.names = F, quote = F, col.names = F)
#     
#     write.table(validation_ids, 
#                 paste0(main_filepath, "clustering/",
#                        p, "_", sx, "/ids_validation.txt"),
#                 sep = "\t", row.names = F, quote = F, col.names = F)
#   })
# })

# Main script ----

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
parser$add_argument("--K", required = TRUE,
                    default = 4,
                    help = "Number of clusters")
parser$add_argument("--L", required = TRUE,
                    default = 2, 
                    help = "Only sample from individuals with at least L measurements")
parser$add_argument("--M", required = TRUE,
                    default = 2,
                    help = "Initialise clusters with diff. at M years post-baseline")
args <- parser$parse_args()

PHENO <- args$phenotype
cat(paste0("PHENO: ", PHENO, "\n"))
SEX_STRATA <- args$ss
cat(paste0("SEX_STRATA: ", SEX_STRATA, "\n"))
K <- as.numeric(args$K)
cat(paste0("number clusters: ", K, "\n"))

main_filepath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/"

# For sampling scheme
NSAMPLES <- 5000 # number of individuals to sample each iteration
S <- 10 # number of times to draw samples from population
L <- as.numeric(args$L) # only sample from individuals with at least L measurements
cat(paste0("minimum # measurements: ", L, "\n"))

# For clustering
M <- args$M
if (M != "random") {
  M <- as.numeric(args$M) # initialise clusters with K-tile diff. at M years post-baseline
}
cat(paste0("K-tile difference at M yrs post-baseline: ", M, "\n"))

resdir <- paste0(main_filepath, "clustering/", PHENO, "_", SEX_STRATA, 
                 "/parameter_selection/")
dir.create(resdir)
plotdir <- paste0(main_filepath, "clustering/", PHENO, "_", SEX_STRATA, 
                 "/parameter_selection/plots/")
dir.create(plotdir)

custom_four_diverge <- c("#D35C79", "#D9AB90", "#9FBCA4", "#009593")
usecolpal <- colorRampPalette(custom_four_diverge)(K)
names(usecolpal) <- paste0("k", 1:K)

# Load data ----

model_dat <- readRDS(paste0(main_filepath, "results/with_rvar_fit_objects_", 
                            PHENO, "_", SEX_STRATA, ".rds"))

covars_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[[PHENO]]
cat(paste0("Number of ids in covariates data: ", nrow(covars_dat), "\n"))

# Only retain 80% of data used for discovery
training_ids <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/", 
                                  PHENO, "_", SEX_STRATA, "/ids_training.txt"),
                           sep = "\t", header = F, stringsAsFactors = F)$V1
training_ids <- as.character(training_ids)

B <- model_dat$B
spline_posteriors <- model_dat$spline_posteriors[training_ids]
model_resid_var <- model_dat$resid_var

# Only retain data with minimum L measurements 
covars_dat <- covars_dat %>% 
  mutate(eid = as.character(eid)) %>%
  filter(eid %in% training_ids & FU_n >= L)
VALID_IDS <- unique(covars_dat$eid)
cat(paste0("\t", "Number of ids in training data: ", length(VALID_IDS), "\n"))

spline_posteriors <- spline_posteriors[VALID_IDS]

# Create B-spline coefficient matrices for distance matrix construction ----

# Create D-matrix for removing intercept
D <- diag(x = 1, nrow = (ncol(B)-1), ncol = (ncol(B)-1))
D <- cbind(rep(-1, (ncol(B)-1)), D)
D <- rbind(rep(0, ncol(B)), D)

# Create mean and S.D. matrix
# Matrix of means (id x basis)
mn_mat <- lapply(spline_posteriors, function (spobj) {
  return (as.data.frame(t(spobj$mu)))
})
mn_mat <- bind_rows(mn_mat)
rownames(mn_mat) <- names(spline_posteriors)
# Remove intercepts
dmn_mat <- t(D %*% t(mn_mat))
rownames(dmn_mat) <- names(spline_posteriors)

# Matrix of standard deviations (id x basis)
sd_mat <- lapply(spline_posteriors, function (spobj) {
  vars <- spobj$Sig
  return (sqrt(diag(vars * model_resid_var)))
})
sd_mat <- as.data.frame(bind_rows(sd_mat))
rownames(sd_mat) <- names(spline_posteriors)

# with intercept removed
dsd_mat <- lapply(spline_posteriors, function (spobj) {
  dvars <- D %*% spobj$Sig %*% t(D)
  return (sqrt(diag(dvars * model_resid_var)))
})
dsd_mat <- t(bind_rows(dsd_mat))
rownames(dsd_mat) <- names(spline_posteriors)

# Function for distance matrix calculation ----

# Given vector of ids, calculate distances between each pair of ids in vector
getCustDistMat <- function (id_list) {
  cat("######## Building distance matrix", "\n")
  sub_mn <- dmn_mat[id_list, ]
  sub_sd <- dsd_mat[id_list, ]
  # Initialise empty distance matrix
  nids <- length(id_list)
  dist_res <- matrix(NA, nids, nids,
                     dimnames = list(id_list, id_list))
  # Loop through each id except the last (as we would have already gone through it)
  for (j in 1:(nids-1)) {
    mn_curr <- sub_mn[j, ]
    mn_diff <- t(apply(sub_mn[j:nids, ], 1, function (x) (x - mn_curr)^2))
    
    v_curr <- sub_sd[j, ]^2
    v_all <- t(apply(sub_sd[j:nids, ], 1, function (x) x^2 + v_curr))
    
    cust_ds <- mn_diff / v_all
    # Replace NaNs that come from 0/0 division with 0s
    cust_ds[is.na(cust_ds)] <- 0
    dist_curr <- sqrt(rowSums(cust_ds))
    dist_res[j, j:nids] <- dist_res[j:nids, j] <- dist_curr
  }
  # Return matrix, replacing the last NA with 0 (as this is the distance of the last id to itself)
  dist_res[nids, nids] <- 0
  return (dist_res)
}

# Functions for clustering ----

# Given a list of ids, calculate medoid (id with min rowsum or colsum 
# in custom distance matrix)
getMedoid <- function (id_list, dmat) {
  sub_mat <- dmat[id_list, ]
  medoid_id <- id_list[which.min(rowSums(sub_mat))]
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
  id_means <- dmn_mat[id_list, ]
  medoid_pos <- dmn_mat[medoid_list, ]
  
  # Euclidean distance from each id to centroid
  closest_medoid <- lapply(1:nrow(id_means), 
                           function (ii) getClosestMedoid(id_means[ii, ], medoid_pos))
  
  res <- data.frame(eid = rownames(id_means),
                    init_clust = unlist(closest_medoid))
  
  return (res)
}

## Calculate cluster medoids for initialisation, given list of individuals 
# and distance matrix for individuals
getInitialCentres <- function (id_list, dmat,
                               myrs = M,
                               ncentres = K) {
  
  cat("######## Initialising cluster centres", "\n")
  
  if (M == "random") {
    
    # Sample initial medoids randomly from list
    init_meds <- sample(id_list, K, replace = F)
    index_ret <- match(init_meds, id_list)
    
    # Assign remaining ids to closest initial centre
    og_assignment <- getClusterBelonging(id_list, init_meds)
    
  } else {
    # Predict full trajectory for all individuals (returns id x time matrix)
    pred_mat <- t(apply(dmn_mat[id_list, ], 1, 
                        function (x) model_dat$B %*% x))
    
    # Get fold-change between 1st and M*365th measurement for each individual
    # and assign to an n-tile
    if (ncol(pred_mat) < myrs*365) stop("Set smaller M")
    diffs_assess <- (pred_mat[, myrs*365] - pred_mat[, 1]) / pred_mat[, 1]
    get_id_q <- data.frame(eid = rownames(pred_mat),
                           fc = diffs_assess) %>%
      mutate(init_clust = ntile(fc, ncentres))
    
    init_meds <- get_id_q %>% 
      group_by(init_clust) %>%
      summarise(gp_medoid = getMedoid(eid, dmat))
    
    # Return indices of medoids and initial cluster assignment
    og_assignment <- get_id_q[, c("eid", "init_clust")]
    index_ret <- match(init_meds$gp_medoid, id_list)
    
  }
  
  return (list(og_assignment = og_assignment,
               medoid_index = index_ret))
}

## Perform PAM clustering (partitioning around medoids) and return clusters

performClustering <- function (dist_mat, ncentres = K, init_centres) {
  clust_res <- pam(dist_mat, k = ncentres, diss = T, 
                   medoids = init_centres, 
                   keep.diss = T, cluster.only = F)
  sil_scores <- silhouette(clust_res)
  
  res <- data.frame(eid = row.names(sil_scores),
                    final_clust = sil_scores[, "cluster"])
  cat("######## Clustering done!", "\n", "\n")
  return (list(clust_assignment = res,
               silhouette_score = mean(sil_scores[, "sil_width"])))
}

# Apply clustering S times ----

itered_clustering <- lapply(1:S, function (si) {
  cat(paste0("Running iteration #", si), "\n")
  # Get ids
  id_list <- sample(VALID_IDS, NSAMPLES, replace = F)
  # Distance matrix for these IDs
  dmat <- getCustDistMat(id_list)
  
  # Initialisation
  init_centres <- getInitialCentres(id_list, dmat,
                                    myrs = M, ncentres = K)
  # Cluster
  clust_res <- performClustering(dist_mat = dmat, 
                                 ncentres = K, 
                                 init_centres = init_centres$medoid_index)
  return (clust_res)
})
names(itered_clustering) <- paste0("iter", 1:S)

# Functions to correspond cluster centroids across iterations ----

# Given a list of grouped ids, get the group centroids and centroid trajectories
# Returns matrix of NCLUST x NDF 
calcGroupCentroids <- function (grouped_id_list) {
  # Get coefficients for the given ids
  sub_mn <- as.data.frame(dmn_mat[grouped_id_list$eid, ])
  sub_mn$eid <- rownames(sub_mn)
  sub_mn$clust <- 
    grouped_id_list$final_clust[match(sub_mn$eid, grouped_id_list$eid)]
  
  # Calculate centroid means
  centroid_means <- sub_mn %>% group_by(clust) %>%
    summarise(across(-eid, mean))
  centroid_means <- centroid_means[, -1]
  
  return (centroid_means)
}

# Given centroid trajectories, reorder from 1 (highest) to NCLUST (lowest)
# Since the trajectories are non-overlapping, just order them by last pred value
reorderCentroids <- function (centroid_traj) {
  new_order <- order(centroid_traj[, ncol(centroid_traj)], decreasing = T)
  map_res <- data.frame(old_clust = 1:nrow(centroid_traj),
                        new_clust = new_order)
  return (map_res)
}

# Plotting functions ----

## Predictions
getPredValuesClusterCentroid <- function (coef_mat) {
  to_calc <- coef_mat[, -1]
  pred_vals <- t(apply(to_calc, 1, 
                       function (x) model_dat$B %*% x))
  # Wrangle into ggplot format
  for_plot <- as.data.frame(pred_vals)
  colnames(for_plot) <- paste0("d", 1:ncol(for_plot))
  for_plot$final_clust <- coef_mat$final_clust
  
  for_plot <- for_plot %>% pivot_longer(cols = -final_clust,
                                        names_to = "t_diff", 
                                        names_prefix = "d", 
                                        values_to = "pred_value") %>%
    mutate(t_diff = as.numeric(t_diff))
  return (for_plot)
}

## Cluster centroid mean and S.D., wrangled into long format for plot
getPlotDatClusterCentroids <- function (grouped_id_list) {
  
  full_dat <- as.data.frame(dmn_mat[grouped_id_list$eid, ])
  full_dat$eid <- rownames(full_dat)
  full_dat$final_clust <- 
    grouped_id_list$final_clust[match(full_dat$eid, grouped_id_list$eid)]
  full_dat$final_clust <- paste0("k", full_dat$final_clust)
  
  centroid_means <- full_dat %>% group_by(final_clust) %>%
    summarise(across(-eid, mean))
  centroid_sds <- full_dat %>% group_by(final_clust) %>%
    summarise(across(-eid, sd))
  
  forloci <- centroid_means[, -1] - 1.96*centroid_sds[, -1]
  forloci$final_clust <- centroid_means$final_clust
  forloci <- forloci[, c("final_clust", paste0("V", 1:(ncol(forloci)-1)))]
  
  forupci <- centroid_means[, -1] + 1.96*centroid_sds[, -1]
  forupci$final_clust <- centroid_means$final_clust
  forupci <- forupci[, c("final_clust", paste0("V", 1:(ncol(forupci)-1)))]
  
  # Predict full trajectory for centroids (returns centroid x time matrix)
  pred_mns <- getPredValuesClusterCentroid(centroid_means) %>%
    rename(pred_mean_value = pred_value)
  pred_locis <- getPredValuesClusterCentroid(forloci) %>%
    rename(pred_loci_value = pred_value)
  pred_upcis <- getPredValuesClusterCentroid(forupci) %>%
    rename(pred_upci_value = pred_value)
  
  # Wrangle into ggplot format
  for_plot <- full_join(pred_mns, pred_locis, by = c("final_clust", "t_diff"))
  for_plot <- full_join(for_plot, pred_upcis, by = c("final_clust", "t_diff"))
  
  return (for_plot)
}

## Confusion matrix plots for initial vs final cluster assignment

getPlotDatConfusionMat <- function (cmat) {
  cmat <- cmat %>% mutate(init_clust = factor(paste0("k", init_clust), 
                                              levels = paste0("k", 1:K)),
                          final_clust = factor(paste0("k", final_clust), 
                                               levels = paste0("k", K:1)))
  res <- ggplot(cmat, aes(x = init_clust, y = final_clust, fill = count)) +
    geom_tile() + 
    geom_text(aes(label = count), size = 2) +
    scale_fill_gradient(low = "white", high = "#009194", guide = "none") +
    labs(x = "Initial", y = "Final") 
  return (res)
}

# Apply to match cluster centroids across iterations ----

cat("Corresponding clusters across iterations", "\n")

combined_iters <- lapply(1:length(itered_clustering), function (si) {
  # Get cluster assignment
  df <- itered_clustering[[paste0("iter", si)]]$clust_assignment
  # Calculate old centroids
  old_centroid_pos <- calcGroupCentroids(df)
  
  # Get old trajectories to reorder
  old_coefs <- t(apply(old_centroid_pos, 1, 
                       function (x) B %*% x))
  # Get new cluster order
  remapping_order <- reorderCentroids(old_coefs)
  # Reassign centroids
  new_centroid_pos <- old_centroid_pos[remapping_order$new_clust, ]
  
  reord <- itered_clustering[[paste0("iter", si)]]$clust_assignment
  reord$final_clust <- 
    remapping_order$old_clust[match(reord$final_clust, remapping_order$new_clust)]
  
  # Pivot longer for tidyverse mean calculations
  for_mean_calc <- as.data.frame(new_centroid_pos)
  colnames(for_mean_calc) <- paste0("b", 1:ncol(for_mean_calc))
  for_mean_calc$clust <- factor(as.character(1:nrow(for_mean_calc)))
  
  for_mean_calc <- for_mean_calc %>% pivot_longer(cols = -clust,
                                                  names_to = "b_coef", 
                                                  names_prefix = "b", 
                                                  values_to = "b_val") %>%
    mutate(b_coef = as.numeric(b_coef),
           iteration = si)
  return (list(reordered = reord,
          coef_dat = for_mean_calc))
})

across_iters_mean_calc <- lapply(combined_iters, function (dlist) {
  return (dlist$coef_dat)
})
across_iters_mean_calc <- bind_rows(across_iters_mean_calc)

# Calculate mean centroid positions across iterations
combined_centroid_mn <- across_iters_mean_calc %>% 
  group_by(clust, b_coef) %>%
  summarise(b_val = mean(b_val))
# Pivot wider back to NCLUST X NDF format
combined_centroid_mn <- combined_centroid_mn %>%
  pivot_wider(id_cols = clust,
              names_from = b_coef, values_from = b_val)

# Calculate mean and S.D. of silhouette score across iterations

combined_silscore <- lapply(1:length(itered_clustering), function (si) { 
  # Get silhouette score
  sscore <- itered_clustering[[paste0("iter", si)]]$silhouette_score
  return (sscore)
})
mean_silscore <- mean(unlist(combined_silscore))
sd_silscore <- sd(unlist(combined_silscore))

silscore <- data.frame(mean = mean_silscore, sd = sd_silscore)

# RETURN: K, L, M, cluster centroids, and silhouette score information ----

to_ret <- list(K = K, L = L, M = M,
               cluster_centroids = combined_centroid_mn,
               silhouette_score = silscore)

saveRDS(to_ret, 
        paste0(resdir, "K", K, "_L", L, "_M", M, ".rds"))

# Plot results ----

## Centroid mean trajectories ---
centroid_traj <- lapply(1:S, function (si) {
  plot_res <- getPlotDatClusterCentroids(combined_iters[[si]]$reordered)
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
  scale_color_manual(values = usecolpal, guide = "none") +
  scale_fill_manual(values = usecolpal, guide = "none") +
  labs(x = "Days from first measurement", 
       y = "Cluster centroid predicted value")

png(paste0(plotdir, "sample_clustering_scheme_", PHENO, "_", SEX_STRATA, 
           "_K", K, "_L", L, "_M", M, 
           ".png"))
print(centroid_traj_plot)
dev.off()

