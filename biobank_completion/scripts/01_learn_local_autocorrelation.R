source("scripts/00_load_data.R")

##################################################################
# Process command line arguments
##################################################################
args_test_presence <- commandArgs(trailingOnly = TRUE)
if (length(args_test_presence) > 0) {
  parser <- ArgumentParser()
  parser$add_argument("--task_id", required = TRUE,
                      help = "Task ID in array")
  parser$add_argument("--n_tasks", required = TRUE,
                      help = "Number of tasks in array")
  args <- parser$parse_args()
  TASK_ID <- as.numeric(args$task_id)
  N_TASKS <- as.numeric(args$n_tasks)
} else {
  TASK_ID <- 1
  N_TASKS <- 1000
}


##############################################################################
# Learn local autocorrelation structure 
##############################################################################
B <- splines::bs(day_unique, df = control$NDF_SPLINE, intercept = TRUE)
ALL_IDS <- names(model_dat)
Sig_smooth <- ar1_covariance(n_time = control$NDF_SPLINE, 
                             rho = control$AR1_RHO, 
                             noise_sd = control$AR1_NOISE_SD, 
                             intercept_sd = control$AR1_INTERCEPT_SD)
precision_smooth <- solve(Sig_smooth)


mean_transform <- matrix(0, N_AGE_BIN, N_DAY)
for (i in 1:N_AGE_BIN) {
  mean_transform[i, day_unique >= age_bin_df$low_day[i] & day_unique <= age_bin_df$upp_day[i]] <- 1
}
mean_transform <- mean_transform / pmax(1, rowSums(mean_transform))


joint_dist_mat <- matrix(0, N_Y_BIN, N_Y_BIN)
joint_dist_mat_multi_list <- list()
for (j in 1:control$N_GRAD_CLASSES) {
  joint_dist_mat_multi_list[[j]] <- rep(list(joint_dist_mat), N_AGE_BIN)
}
force_autocorr_calc <- FALSE
subj_split <- subj_split_fn(all_subj, N_TASKS = N_TASKS)  
subj_do <- subj_split[[TASK_ID]]
for (subj_sim in subj_do) {print(subj_sim)
  index_subj <- match(subj_sim, all_subj)
  spline_posterior_curr <- lapply(model_dat[subj_sim], function (id_df) {
    y <- id_df$value
    X <- matrix(0, nrow = length(day_unique), ncol = nrow(id_df))
    rownames(X) <- day_unique
    for (j in 1:nrow(id_df)) {
      fill_val <- as.numeric(id_df[j, "age_days"])
      X[as.character(fill_val), j] <- 1
    }
    res <- fit_subj_posterior_under_simple_prior(Z = t(B) %*% X,
                                                 y = y,
                                                 precision_smooth = precision_smooth)
    return (res)
  })
  
  post_spline_coefs <- t(MASS::mvrnorm(n = control$N_MC_SAMPLES, mu = spline_posterior_curr[[subj_sim]]$mu, 
                                       Sigma = spline_posterior_curr[[subj_sim]]$Sig * control$RESID_NOISE_VAR)) 
  sims <- mean_transform %*% B %*% post_spline_coefs
  y_bin_sims <- sims
  y_bin_sims[] <- findInterval(x = sims, vec = y_bin_df$l) 
  
  curr_subj_time_range <- subj_time_ranges[match(subj_sim, subj_time_ranges$eid), ]
  curr_subj_bins_update_logical <- sapply(1:N_AGE_BIN, function(i) curr_subj_time_range$lower <= age_bin_df$upp_day[max(c(1, i - 2))] & 
                                            curr_subj_time_range$upper >= age_bin_df$low_day[min(c(N_AGE_BIN, i))] & i >= 3)
  for (age_bin_curr in which(curr_subj_bins_update_logical)) {
    joint_dist <- t(y_bin_sims[age_bin_curr + (-2):0, ])
    for (j in 1:nrow(joint_dist)) {
      coords_curr <- joint_dist[j, ]
      class_vec <- (1:control$N_GRAD_CLASSES) - mean(1:control$N_GRAD_CLASSES)
      grad_curr <- coords_curr[2] - coords_curr[1]
      grad_curr <- max(c(min(class_vec), grad_curr))
      grad_curr <- min(c(max(class_vec), grad_curr))
      list_to_add_to <- grad_curr + mean(1:control$N_GRAD_CLASSES)
      joint_dist_mat_multi_list[[list_to_add_to]][[age_bin_curr]][joint_dist[j, 2], joint_dist[j, 3]] <- 
        joint_dist_mat_multi_list[[list_to_add_to]][[age_bin_curr]][joint_dist[j, 2], joint_dist[j, 3]] + 1
    }
  }
}
dir.create("output/joint_dist_mat", showWarnings = FALSE)
saveRDS(joint_dist_mat_multi_list, file = paste0("output/joint_dist_mat/joint_dist_mat_multi_list_task_", TASK_ID, ".RDS"))







