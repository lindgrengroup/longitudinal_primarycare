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

RANDOM_SEED <- 160522
set.seed(RANDOM_SEED)

parser <- ArgumentParser()
parser$add_argument("--phenotype", required = TRUE,
                    help = "Phenotype to model")
parser$add_argument("--ss", required = TRUE,
                    help = "Sex strata")
args <- parser$parse_args()

PHENO <- args$phenotype
cat(paste0("PHENO: ", PHENO, "\n"))
SEX_STRATA <- args$ss
cat(paste0("SEX_STRATA: ", SEX_STRATA, "\n"))
K <- 4
cat(paste0("number clusters: ", K, "\n"))

main_filepath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/"

# For sampling scheme
NSAMPLES <- 5000 # number of individuals to sample each iteration
S <- 10 # number of times to draw samples from population
L <- 2 # only sample from individuals with at least L measurements
cat(paste0("minimum # measurements: ", L, "\n"))

# For clustering
M <- 2
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

HOLDOUT_IDS <- sample(training_ids, NSAMPLES, replace = F)
FOR_CLUSTERING_IDS <- training_ids[!training_ids %in% HOLDOUT_IDS]

B <- model_dat$B
spline_posteriors <- model_dat$spline_posteriors[training_ids]
model_resid_var <- model_dat$resid_var

# Only retain data with minimum L measurements 
covars_dat <- covars_dat %>% 
  mutate(eid = as.character(eid)) %>%
  filter(eid %in% training_ids & FU_n >= L)
VALID_IDS <- unique(covars_dat$eid)

HOLDOUT_IDS <- HOLDOUT_IDS[HOLDOUT_IDS %in% VALID_IDS]
FOR_CLUSTERING_IDS <- FOR_CLUSTERING_IDS[FOR_CLUSTERING_IDS %in% VALID_IDS]

cat(paste0("\t", "Number of ids in training data: ", length(FOR_CLUSTERING_IDS), "\n"))

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
  id_list <- sample(FOR_CLUSTERING_IDS, NSAMPLES, replace = F)
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

# Apply to match cluster centroids across iterations ----

cat("Corresponding clusters across iterations", "\n")

remapped_centroids <- lapply(1:S, function (si) {
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
  return (new_centroid_pos)
})

# Assign holdout ids to clusters from each random split ----

holdout_mns <- dmn_mat[HOLDOUT_IDS, ]

holdout_assignments <- lapply(1:S, function (si) {
  medoid_pos <- remapped_centroids[[si]]
  
  # Euclidean distance from each id to centroid
  closest_medoid <- lapply(1:nrow(holdout_mns), 
                           function (ii) getClosestMedoid(holdout_mns[ii, ], medoid_pos))
  
  ret_df <- data.frame(eid = HOLDOUT_IDS,
                       clust = paste0("k", closest_medoid))
  colnames(ret_df) <- c("eid", paste0("clust_iter", si))
  return (ret_df)
})
holdout_assignments <- holdout_assignments %>%
  reduce(full_join, by = "eid")

write.table(holdout_assignments,
            paste0(resdir, "holdout_cluster_assignments_", PHENO, "_", SEX_STRATA, 
                   "_K", K, "_L", L, "_M", M, 
                   ".txt"),
            sep = "\t", row.names = F, quote = F)

# Get modal cluster for each ID and number of times cluster is assigned to mode

getMode <- function (x) {
  ux <- unique(x)
  modal_clust <- ux[which.max(tabulate(match(x, ux)))]
  return (modal_clust)
}

modal_clust_per_id <- apply(holdout_assignments[,-1], 1, 
                            function (x) getMode(x))
freq_mode <- lapply(1:nrow(holdout_assignments), function (i) {
  x <- holdout_assignments[i, -1]
  return (sum(x == modal_clust_per_id[i]))
})
freq_mode <- unlist(freq_mode)

holdout_assignments$modal_clust <- modal_clust_per_id
holdout_assignments$freq_mode <- freq_mode

# Plot results ----

resplot <- ggplot(holdout_assignments, aes(x = freq_mode,
                                           fill = modal_clust, colour = modal_clust)) +
  facet_wrap(~modal_clust, nrow = 2) +
  geom_histogram(binwidth = 1, position = "identity") +
  scale_color_manual(values = usecolpal, guide = "none") +
  scale_fill_manual(values = usecolpal, guide = "none") +
  labs(x = "Number of times hold-out id assigned to its modal cluster",
       y = "Number of ids")

png(paste0(plotdir, "holdout_cluster_assignments_", PHENO, "_", SEX_STRATA, 
           "_K", K, "_L", L, "_M", M, 
           ".png"))
print(resplot)
dev.off()

