# Author: Samvida S. Venkatesh
# Date: 19/08/2022

library(argparse)
library(MASS)
library(tidyverse)

# FROM PARAMETER SELECTION PLOTS, WE SEE THAT 
# K = 4 (number of clusters), 
# L = 2 (minimum number of follow-up measures),
# M = 2 (initialisation at obesity-change 2 years post-baseline)
# yield the most dense and separable clusters
K_chosen <- 4
L_chosen <- 2
M_chosen <- 2

RANDOM_SEED <- 160522
set.seed(RANDOM_SEED)

# Get arguments ----

parser <- ArgumentParser()
parser$add_argument("--phenotype", required = TRUE,
                    help = "Phenotype to model")
parser$add_argument("--ss", required = TRUE,
                    help = "Sex strata")
parser$add_argument("--nboots", default = 100,
                    required = TRUE,
                    help = "Number of bootstrap samples to draw per individual")

args <- parser$parse_args()

PHENO <- args$phenotype
SEX_STRATA <- args$ss
NBOOTS <- as.numeric(args$nboots)

resdir <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/clustering/", 
                 PHENO, "_", SEX_STRATA)

# Load data ----

model_dat <- readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/results/with_rvar_fit_objects_", 
                            PHENO, "_", SEX_STRATA, ".rds"))

clust_centres <- readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/clustering/", 
                                PHENO, "_", SEX_STRATA, "/parameter_selection/K",
                                K_chosen, "_L", L_chosen, "_M", M_chosen, ".rds"))
CLUST_NAMES <- paste0("k", 1:K_chosen)

# Calculate mean and S.D. matrix of coefficients needed for distance calculations ----

B <- model_dat$B
spline_posteriors <- model_dat$spline_posteriors
model_resid_var <- model_dat$resid_var

# Create D-matrix for removing intercept
D <- diag(x = 1, nrow = (ncol(B)-1), ncol = (ncol(B)-1))
D <- cbind(rep(-1, (ncol(B)-1)), D)
D <- rbind(rep(0, ncol(B)), D)

# Matrix of means (id x basis)
mn_mat <- lapply(spline_posteriors, function (spobj) {
  return (as.data.frame(t(spobj$mu)))
})
mn_mat <- bind_rows(mn_mat)
rownames(mn_mat) <- names(spline_posteriors)
# Remove intercepts
dmn_mat <- t(D %*% t(mn_mat))
rownames(dmn_mat) <- names(spline_posteriors)

# Functions to generate bootstrapped position vectors for an individual ----

createCoefSamples <- function (id, nboots = 100) {
  mu_vec <- dmn_mat[id, ]
  cov_mat <- D %*% spline_posteriors[[id]]$Sig %*% t(D)
  res <- mvrnorm(n = nboots,
                 mu = mu_vec,
                 Sigma = cov_mat)
  return (res)
}

# Functions to assign position vectors to their closest cluster ----

# Given a vector and a matrix, get row index in matrix closest to the vector
getClosestMedoid <- function (x1, mat_xs) {
  x1 <- unlist(x1)
  all_dists <- t(apply(mat_xs, 1, function (x) x - x1))^2
  return (which.min(rowSums(all_dists)))
}

# Given a matrix of id-positions and centroid positions, assign each position
# to closest medoid
# Return probability of belonging to each cluster
getClusterProbs <- function (posn_mat, centroid_mat) {
  # Based on Euclidean distance from each id to centroid
  closest_medoid <- sapply(1:nrow(posn_mat), 
                           function (ii) getClosestMedoid(posn_mat[ii, ], centroid_mat))
  count_probs <- as.data.frame(table(closest_medoid)) %>%
    mutate(closest_medoid = paste0("k", closest_medoid),
           prob = Freq / nrow(posn_mat))  
  return (count_probs[, c("closest_medoid", "prob")])
}

# Get cluster probabilities ----

ALL_IDS <- names(spline_posteriors)
# First column of centroid matrix is cluster name
centroid_mat <- clust_centres$cluster_centroids[, -1]

soft_clust_res <- lapply(ALL_IDS, function (id) {
  bsamples <- createCoefSamples(id, nboots = NBOOTS)
  cluster_probs <- getClusterProbs(bsamples, centroid_mat) 
  cluster_probs$eid <- id
  return (cluster_probs)
})
soft_clust_res <- bind_rows(soft_clust_res) %>%
  pivot_wider(id_cols = eid,
              names_from = closest_medoid,
              values_from = prob, values_fill = 0)
soft_clust_res <- soft_clust_res[, c("eid", CLUST_NAMES)]

write.table(soft_clust_res,
            paste0(resdir, "/soft_clustering_probs_", PHENO, "_", SEX_STRATA, 
                   ".txt"), sep = "\t", row.names = F, quote = F)
