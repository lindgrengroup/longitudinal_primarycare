# Author: Samvida S. Venkatesh
# Date: 24/05/2022

library(splines)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
theme_set(theme_bw())

RANDOM_SEED <- 160522
set.seed(RANDOM_SEED)

# Get arguments ----

parser <- ArgumentParser()
parser$add_argument("--phenotype", required = TRUE,
                    help = "Phenotype to model")
parser$add_argument("--sex_strata", required = TRUE,
                    help = "Sex strata")

args <- parser$parse_args()

PHENO <- args$phenotype
SEX_STRATA <- args$sex_strata

plotdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/plots/sample_clustering_scheme/"
resdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/"

# Load data ----

model_dat <- readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/results/fit_objects_", 
                            PHENO, "_", SEX_STRATA, ".rds"))

model_dat$resid_var <- 5
saveRDS(model_dat, 
        paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/results/fit_objects_", 
               PHENO, "_", SEX_STRATA, ".rds"))

L <- c(2, 5, 10)
M <- c("random", 1, 2, 5, 10)

# Results are across multiple iterations so access with res[[lmin]][[myrs]][[iter]]
clust_res <- lapply(L, function (lmin) {
  res_list <- lapply(M, function (myrs) {
    readRDS(paste0(resdir, "sample_clustering_scheme_", PHENO, "_", SEX_STRATA, 
           "_L", lmin, "_M", myrs, 
           ".rds"))
  })
  names(res_list) <- paste0("M", M)
  return (res_list)
})
names(clust_res) <- paste0("L", L)

# Calculate mean and S.D. matrix of coefficients needed for distance calculations ----

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

# Functions to combine cluster centroids across iterations ----

# Given a list of grouped ids, get the group centroids and centroid trajectories
# Returns matrix of NCLUST x NDF 
calcCentroids <- function (grouped_id_list) {
  # Get coefficients for the given ids
  sub_mn <- mn_mat[grouped_id_list$eid, ]
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

# Apply to correspond cluster centroids and get mean of mean centroid ----

# Return NCLUST x NDF matrix with final centroid positions
final_centroid_pos <- lapply(L, function (lmin) {
  res_list <- lapply(M, function (myrs) {
    iter_list <- clust_res[[paste0("L", lmin)]][[paste0("M", myrs)]]
    
    iter_res <- lapply(1:length(iter_list), function (si) {
      # Get cluster assignment
      df <- iter_list[[paste0("iter", si)]]
      # Calculate old centroids
      old_centroid_pos <- calcCentroids(df)
      
      # Get old trajectories to reorder
      old_coefs <- t(apply(old_centroid_pos, 1, 
                            function (x) B %*% x))
      # Get new cluster order
      remapping_order <- reorderCentroids(old_coefs)
      # Reassign centroids
      new_centroid_pos <- old_centroid_pos[remapping_order$new_clust, ]
      
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
      return (for_mean_calc)
    })
    iter_res <- bind_rows(iter_res)
    
    # Calculate mean centroid positions across iterations
    # Also output S.D.s for plot
    combined_centroid_mn <- iter_res %>% 
      group_by(clust, b_coef) %>%
      summarise(b_val = mean(b_val))
    # Pivot wider back to NCLUST X NDF format
    combined_centroid_mn <- combined_centroid_mn %>%
      pivot_wider(id_cols = clust,
                  names_from = b_coef, values_from = b_val)
    
    combined_centroid_sd <- iter_res %>% 
      group_by(clust, b_coef) %>%
      summarise(b_val = sd(b_val))
    # Pivot wider back to NCLUST X NDF format
    combined_centroid_sd <- combined_centroid_sd %>%
      pivot_wider(id_cols = clust,
                  names_from = b_coef, values_from = b_val)
    
    return (list(mean_mat = combined_centroid_mn[, -1],
                 sd_mat = combined_centroid_sd[, -1]))
  })
  names(res_list) <- paste0("M", M)
  return (res_list)
})
names(final_centroid_pos) <- paste0("L", L)

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

## Cluster centroid mean and S.D., wrangled into long format for plot
getPlotDatClusterCentroids <- function (centroid_pos_mats) {
  
  # Predict full trajectory for centroids (returns centroid x time matrix)
  pred_mns <- getPredValuesClusterCentroid(centroid_pos_mats$mean_mat) %>%
    rename(pred_mean_value = pred_value)
  pred_locis <- getPredValuesClusterCentroid(centroid_pos_mats$mean_mat - 1.96*centroid_pos_mats$sd_mat) %>%
    rename(pred_loci_value = pred_value)
  pred_upcis <- getPredValuesClusterCentroid(centroid_pos_mats$mean_mat + 1.96*centroid_pos_mats$sd_mat) %>%
    rename(pred_upci_value = pred_value)
  
  # Wrangle into ggplot format
  for_plot <- full_join(pred_mns, pred_locis, by = c("clust", "t_diff"))
  for_plot <- full_join(for_plot, pred_upcis, by = c("clust", "t_diff"))
  
  return (for_plot)
}

## Apply ----

sanity_check_cluster_centres <- lapply(L, function (lmin) {
  mlist <- lapply(M, function (myrs) {
    cluster_centres <- final_centroid_pos[[paste0("L", lmin)]][[paste0("M", myrs)]]
    plot_dat <- getPlotDatClusterCentroids(cluster_centres) %>%
      mutate(mlname = paste0("L", lmin, "_M", myrs))
    return (plot_dat)
  })
  mlist <- bind_rows(mlist)
  # Plot a column for all M-values for this L-value
  mcol_plot <- ggplot(mlist, aes(x = t_diff, y = pred_mean_value,
                                   col = clust, fill = clust)) +
    facet_wrap(~mlname, ncol = 1) +
    geom_line() +
    geom_ribbon(aes(ymin = pred_loci_value, ymax = pred_upci_value),
                alpha = 0.1, linetype = 0) +
    scale_color_brewer(palette = "Set1", guide = "none") +
    scale_fill_brewer(palette = "Set1", guide = "none") +
    labs(x = "Days from first measurement", 
         y = "Cluster centroid predicted value")
  return (mcol_plot)
})

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

# Go through all possible initialisations and see how clusters are assigned ----

ALL_IDS <- rownames(mn_mat)

all_ids_cluster_assignment <- lapply(L, function (lmin) {
  res_list <- lapply(M, function (myrs) {
    print(paste0("Assigning clusters for Lmin = ", lmin, 
               "; Myrs = ", myrs, "\n"))
    return (getClusterBelonging(ALL_IDS, 
                                final_centroid_pos[[paste0("L", lmin)]][[paste0("M", myrs)]]$mean_mat))
    
  })
  names(res_list) <- paste0("M", M)
  return (res_list)
})
names(all_ids_cluster_assignment) <- paste0("L", L)

saveRDS(all_ids_cluster_assignment,
        paste0(resdir, "all_assigned_clusters_", PHENO, "_", SEX_STRATA, 
               ".rds"))

## Function to plot confusion matrix 
plotConfusionMat <- function (cmat, c1_label, c2_label) {
  nclust <- length(unique(cmat$c1))
  cmat <- cmat %>% mutate(c1 = factor(c1, levels = 1:nclust),
                          c2 = factor(c2, levels = nclust:1))
  res <- ggplot(cmat, aes(x = c1, y = c2, fill = count)) +
    geom_tile() + 
    geom_text(aes(label = count), size = 2) +
    scale_fill_gradient(low = "white", high = "#009194", guide = "none") +
    labs(x = c1_label, y = c2_label) 
  return (res)
}

## Apply plotting ----
mlcombos <- paste0("L", rep(L, each = 5), "_M", rep(M, times = 3))

confusion_plots <- lapply(mlcombos, function (ml) {
  c1 <- ml
  get_l_c1 <- gsub("L", "", gsub("_M.*", "", c1))
  get_m_c1 <- gsub(".*_M", "", c1)
  ml_other <- mlcombos[mlcombos != c1]
  
  plot_list <- lapply(ml_other, function (c2) {
    get_l_c2 <- gsub("L", "", gsub("_M.*", "", c2))
    get_m_c2 <- gsub(".*_M", "", c2)
    
    cmat_c1 <- all_ids_cluster_assignment[[paste0("L", get_l_c1)]][[paste0("M", get_m_c1)]] %>%
      rename(c1 = clust)
    cmat_c2 <- all_ids_cluster_assignment[[paste0("L", get_l_c2)]][[paste0("M", get_m_c2)]] %>%
      rename(c2 = clust)
    
    cmat <- full_join(cmat_c1, cmat_c2, by = "eid")
    cmat <- as.data.frame(table(cmat$c1, cmat$c2))
    colnames(cmat) <- c("c1", "c2", "count")
    return (plotConfusionMat(cmat, c1_label = c1, c2_label = c2))
  })
  return (ggarrange(plotlist = plot_list, 
                    ncol = 3, nrow = 5))
})

pdf(paste0(plotdir, "combined_iterations_clusters_", PHENO, "_", SEX_STRATA, 
           ".pdf"), onefile = T)
# Arrange cluster centroid plots
print(ggarrange(plotlist = sanity_check_cluster_centres, ncol = 3, nrow = 1))
# Print all the confusion matrices
print(confusion_plots)
dev.off()

# Choose which set of L and M to save cluster assignments for ----

to_write <- all_ids_cluster_assignment[["L2"]][["M5"]]
write.table(to_write, paste0(resdir, "assigned_clusters_", PHENO, "_", SEX_STRATA, 
                             ".txt"),
            sep = "\t", quote = F, row.names = F)
