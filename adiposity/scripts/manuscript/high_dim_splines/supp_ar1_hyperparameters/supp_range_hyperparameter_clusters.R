# Author: Samvida S. Venkatesh
# Adapted from: George Nicholson
# Date: 16/05/22

# Main script ----

library(argparse)
library(splines)
library(cluster)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
theme_set(theme_bw())

mainpath <- "" # REDACTED
gen_resources_path <- "" # REDACTED
hidim_mods_path <- "" # REDACTED

RANDOM_SEED <- 160522
set.seed(RANDOM_SEED)

parser <- ArgumentParser()
parser$add_argument("--phenotype", required = TRUE,
                    help = "Phenotype to model")
parser$add_argument("--ss", required = TRUE,
                    help = "Sex strata")
args <- parser$parse_args()

K <- 4
L <- 2
M <- 2

PHENO <- args$phenotype
cat(paste0("PHENO: ", PHENO, "\n"))
SEX_STRATA <- args$ss
cat(paste0("SEX_STRATA: ", SEX_STRATA, "\n"))
cat(paste0("number clusters: ", K, "\n"))
cat(paste0("minimum # measurements: ", L, "\n"))
cat(paste0("K-tile difference at M yrs post-baseline: ", M, "\n"))

resdir <- paste0(hidim_mods_path, "/ar1_parameter_selection/clustering/")
dir.create(resdir)

custom_four_diverge <- c("#D35C79", "#D9AB90", "#9FBCA4", "#009593")
usecolpal <- colorRampPalette(custom_four_diverge)(K)
names(usecolpal) <- paste0("k", 1:K)

# Load data ----

model_dat <- lapply(c(1:3), function (set_i) {
  res <- readRDS(paste0(hidim_mods_path, "/ar1_parameter_selection/fit_objects_", PHENO, "_", SEX_STRATA, 
                 "_parameter_set_", set_i, ".rds"))
  res$spline_posteriors <- res$spline_posteriors[[set_i]]
  saveRDS(res, paste0(hidim_mods_path, "/ar1_parameter_selection/fit_objects_", PHENO, "_", SEX_STRATA, 
                      "_parameter_set_", set_i, ".rds"))
  return (res)
})

# Create B-spline coefficient matrices for distance matrix construction ----

B <- model_dat[[1]]$B
D <- diag(x = 1, nrow = (ncol(B)-1), ncol = (ncol(B)-1))
D <- cbind(rep(-1, (ncol(B)-1)), D)
D <- rbind(rep(0, ncol(B)), D)

dmn_mats <- lapply(c(1:3), function (set_i) {
  spline_posteriors <- model_dat[[set_i]]$spline_posteriors
  model_resid_var <- model_dat[[set_i]]$resid_var
  
  # Matrix of means (id x basis)
  mn_mat <- lapply(spline_posteriors, function (spobj) {
    return (as.data.frame(t(spobj$mu)))
  })
  mn_mat <- bind_rows(mn_mat)
  rownames(mn_mat) <- names(spline_posteriors)
  # Remove intercepts
  dmn_mat <- t(D %*% t(mn_mat))
  rownames(dmn_mat) <- names(spline_posteriors)
  return (dmn_mat)
})

dsd_mats <- lapply(c(1:3), function (set_i) {
  spline_posteriors <- model_dat[[set_i]]$spline_posteriors
  model_resid_var <- model_dat[[set_i]]$resid_var
  
  dsd_mat <- lapply(spline_posteriors, function (spobj) {
    dvars <- D %*% spobj$Sig %*% t(D)
    return (sqrt(diag(dvars * model_resid_var)))
  })
  dsd_mat <- t(bind_rows(dsd_mat))
  rownames(dsd_mat) <- names(spline_posteriors)
  return (dsd_mat)
})

# Function for distance matrix calculation ----

# Given vector of ids, calculate distances between each pair of ids in vector
getCustDistMat <- function (dmn_mat, dsd_mat) {
  cat("######## Building distance matrix", "\n")
  # Initialise empty distance matrix
  nids <- nrow(dmn_mat)
  dist_res <- matrix(NA, nids, nids,
                     dimnames = list(rownames(dmn_mat), rownames(dmn_mat)))
  # Loop through each id except the last (as we would have already gone through it)
  for (j in 1:(nids-1)) {
    mn_curr <- dmn_mat[j, ]
    mn_diff <- t(apply(dmn_mat[j:nids, ], 1, function (x) (x - mn_curr)^2))
    
    v_curr <- dsd_mat[j, ]^2
    v_all <- t(apply(dsd_mat[j:nids, ], 1, function (x) x^2 + v_curr))
    
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
getMedoid <- function (id_list, dist_mat) {
  sub_mat <- dist_mat[id_list, ]
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
getClusterBelonging <- function (id_list, medoid_list, dmn_mat) {
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
getInitialCentres <- function (dmn_mat, dist_mat,
                               myrs = M,
                               ncentres = K) {
  id_list <- rownames(dmn_mat)
  
  cat("######## Initialising cluster centres", "\n")
  
  if (M == "random") {
    
    # Sample initial medoids randomly from list
    init_meds <- sample(id_list, K, replace = F)
    index_ret <- match(init_meds, id_list)
    
    # Assign remaining ids to closest initial centre
    og_assignment <- getClusterBelonging(id_list, init_meds, dmn_mat)
    
  } else {
    # Predict full trajectory for all individuals (returns id x time matrix)
    pred_mat <- t(apply(dmn_mat[id_list, ], 1, 
                        function (x) B %*% x))
    
    # Get fold-change between 1st and M*365th measurement for each individual
    # and assign to an n-tile
    if (ncol(pred_mat) < myrs*365) stop("Set smaller M")
    diffs_assess <- (pred_mat[, myrs*365] - pred_mat[, 1]) / pred_mat[, 1]
    get_id_q <- data.frame(eid = rownames(pred_mat),
                           fc = diffs_assess) %>%
      mutate(init_clust = ntile(fc, ncentres))
    
    init_meds <- get_id_q %>% 
      group_by(init_clust) %>%
      summarise(gp_medoid = getMedoid(eid, dist_mat))
    
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
  return (res)
}

# Apply clustering ----

per_parameter_clustering <- lapply(c(1:3), function (set_i) {
  cat(paste0("Running parameter set #", set_i), "\n")
  
  # Build distance matrix
  dist_mat <- getCustDistMat(dmn_mats[[set_i]], dsd_mats[[set_i]])
  
  # Initialisation
  init_centres <- getInitialCentres(dmn_mat = dmn_mats[[set_i]], 
                                    dist_mat = dist_mat,
                                    myrs = M, ncentres = K)
  # Cluster
  clust_res <- performClustering(dist_mat = dist_mat, 
                                 ncentres = K, 
                                 init_centres = init_centres$medoid_index)
  return (clust_res)
})

# Order centroids from 1 (highest) to 4 (lowest) ----

# Given a list of grouped ids, get the group centroids and centroid trajectories
# Returns matrix of NCLUST x NDF 
calcGroupCentroids <- function (grouped_id_list, dmn_mat) {
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

final_clusts <- lapply(c(1:3), function (set_i) {
  # Get cluster assignment
  df <- per_parameter_clustering[[set_i]]
  # Calculate old centroids
  old_centroid_pos <- calcGroupCentroids(df, dmn_mat = dmn_mats[[set_i]])
  
  # Get old trajectories to reorder
  old_coefs <- t(apply(old_centroid_pos, 1, 
                       function (x) B %*% x))
  # Get new cluster order
  remapping_order <- reorderCentroids(old_coefs)
  # Reassign centroids
  new_centroid_pos <- old_centroid_pos[remapping_order$new_clust, ]
  
  reord <- per_parameter_clustering[[set_i]]
  reord$final_clust <- 
    remapping_order$old_clust[match(reord$final_clust, remapping_order$new_clust)]
  colnames(reord) <- c("eid", paste0("parameter_set_", set_i))
  
  return (reord)
})

# Plot confusion matrices ----

clust_assts <- final_clusts %>% reduce(full_join, by = "eid")

set1_v2 <- as.data.frame(table(clust_assts$parameter_set_1, clust_assts$parameter_set_2)) %>%
  mutate(set1 = factor(Var1, 
                       levels = paste0("k", 1:K)),
         set2 = factor(Var2, 
                       levels = paste0("k", K:1)))
colnames(set1_v2) <- c("set1", "set2", "count")

res <- ggplot(set1_v2, aes(x = set2, y = set1, fill = count)) +
  geom_tile() + 
  geom_text(aes(label = count), size = 2) +
  scale_fill_gradient(low = "white", high = "#009194", guide = "none") +
  labs(x = "Parameter set 2", y = "Parameter set 1") +
  theme(axis.text = element_text(size = 7))

png(paste0(resdir, "clust_assignment_check_set1_vs_set2_", 
           PHENO, "_", SEX_STRATA, ".png"), res = 300,
    height = 5, width = 5, units = "cm")
print(res)
dev.off()

set3_v2 <- as.data.frame(table(clust_assts$parameter_set_3, clust_assts$parameter_set_2)) %>%
  mutate(set3 = factor(Var1, 
                       levels = paste0("k", 1:K)),
         set2 = factor(Var2, 
                       levels = paste0("k", K:1)))
colnames(set3_v2) <- c("set3", "set2", "count")

res <- ggplot(set3_v2, aes(x = set2, y = set3, fill = count)) +
  geom_tile() + 
  geom_text(aes(label = count), size = 2) +
  scale_fill_gradient(low = "white", high = "#009194", guide = "none") +
  labs(x = "Parameter set 2", y = "Parameter set 3") +
  theme(axis.text = element_text(size = 7))

png(paste0(resdir, "clust_assignment_check_set3_vs_set2_", 
           PHENO, "_", SEX_STRATA, ".png"), res = 300,
    height = 5, width = 5, units = "cm")
print(res)
dev.off()

