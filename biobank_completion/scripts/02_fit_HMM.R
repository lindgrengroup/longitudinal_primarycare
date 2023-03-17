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


trans_mats <- readRDS(file = "output/trans_mats.RDS")


###########################################
# Masking data parameters
###########################################
mask_data <- FALSE
window_widths <- sapply(model_dat, function(d) diff(range(d$age_days)))
# summary(window_widths / 365.25)
window_means <- sapply(model_dat, function(d) mean(d$age_days))
# mean(window_means / 365.25)
AGE_MN_TRAIN <- 57
AGE_RADIUS_TRAIN <- 6
N_STATE <- N_Y_BIN * control$N_GRAD_CLASSES


##############################################################################
# Fit HMM one subj at a time
##############################################################################
subj_split <- subj_split_fn(all_subj, N_TASKS = N_TASKS)  
subj_do <- subj_split[[TASK_ID]]
list_out <- list()
for (subj_test in subj_do) {
  if (control$do_plots) {
    pdf(file = paste0("plots/prelim_lines_set", TASK_ID, ".pdf"), width = 9, height = 9)
  }
  # subj_test <- sample(all_subj, 1)
  dat_all <- as.data.frame(model_dat[[subj_test]])
  dat_all$age_yrs <- dat_all$age_days / 365.25
  dat_all$age_yrs_rnd <- round(dat_all$age_yrs)
  if (mask_data) {
    min_age_train <- (AGE_MN_TRAIN - AGE_RADIUS_TRAIN) * 365.25
    max_age_train <- (AGE_MN_TRAIN + AGE_RADIUS_TRAIN) * 365.25
    dat_train <- dat_all[dat_all$age_days >= min_age_train & dat_all$age_days <= max_age_train, ]
    dat_test <- dat_all[dat_all$age_days < min_age_train | dat_all$age_days > max_age_train, ]
  } else {
    dat_train <- dat_all
  }    
  
  ages_obs <- unique(dat_train$age_yrs_rnd)
  dat_train_age_averaged <- dat_train[match(ages_obs,  dat_train$age_yrs_rnd), ]
  dat_train_age_averaged$y_mn <- sapply(ages_obs, function(x) mean(dat_train$value[dat_train$age_yrs_rnd == x]))
  dat_train_age_averaged$y_n <- sapply(ages_obs, function(x) sum(dat_train$age_yrs_rnd == x))
  
  if (NROW(dat_train) > 0) {
    
    ##############################################################################
    # Calculate likelihood matrix
    ##############################################################################
    log_lik_mat <- matrix(0, N_AGE_BIN, N_Y_BIN)
    for (ind in 1:nrow(dat_train)) {
      age_curr <- dat_train$age_days[ind]
      bin_num_curr <- which(age_bin_df$low_day <= age_curr & age_bin_df$upp_day >= age_curr )
      age_l <- age_bin_df$low_day[bin_num_curr]
      age_u <- age_bin_df$upp_day[bin_num_curr]
      for (j in 1:N_Y_BIN) {
        eval_pts <- seq(y_bin_df[j, "l"], y_bin_df[j, "u"], length.out = control$N_LIK_EVAL_PTS)
        log_lik_vec <- dnorm(x = dat_train$value[ind], mean = eval_pts, sd = control$HMM_RESID_NOISE_VAR, log = TRUE)
        log_integrated_lik <- matrixStats::logSumExp(log_lik_vec) - log(control$N_LIK_EVAL_PTS)
        log_lik_mat[bin_num_curr, j] <- log_lik_mat[bin_num_curr, j] + log_integrated_lik
      }
    }
    lik <- exp(log_lik_mat - apply(log_lik_mat, 1, max))
    ##############################################################################
    # Forwards pass
    ##############################################################################
    parl.zero <- parl.inf <- list()
    for (j in 1:N_AGE_BIN) {
      parl.zero[[j]] <- matrix(0, N_STATE, N_STATE)
      parl.inf[[j]] <- matrix(-Inf, N_STATE, N_STATE)
    }
    parl <- parl.zero
    parll <- parl.inf
    for(j in 1:N_AGE_BIN){
      if(j == 1){
        for (k in 1:N_STATE) {
          parl[[j]][k, ] <- lik[j, ]
        }
        parl[[j]] <- parl[[j]] / sum(parl[[j]])
      }
      if (j > 1) {
        fwd_probs <- colSums(parl[[j - 1]])
        updatec <- which(fwd_probs != 0)
        parll[[j]][updatec, ] <- sweep(log(fwd_probs[updatec]) + log(trans_mats[[j - 1]][updatec, , drop = F]), 2, log_lik_mat[j, ], '+')
        parl[[j]][updatec, ] <- exp(parll[[j]][updatec, , drop = F] - max(parll[[j]][updatec, , drop = F]))
        parl[[j]][updatec, ] <- parl[[j]][updatec, , drop = F] / sum(parl[[j]][updatec, , drop = F])
      }
    }
    
    ##############################################################################
    # Multiple random backwards passes
    ##############################################################################
    out_mat <- NULL
    prob_mat <- 0
    for (it in 1:control$N_ITS) {
      out <- rep(NA, N_AGE_BIN)
      for (j in N_AGE_BIN:1) {
        if (j == N_AGE_BIN) {
          multinom_sam_probs <- colSums(parl[[j]])
        }
        if (j < N_AGE_BIN) {
          multinom_sam_probs <- parl[[j + 1]][, out[j + 1]]
        }
        multinom_sam_probs <- multinom_sam_probs / sum(multinom_sam_probs)
        out[j] <- which(rmultinom(1, 1, multinom_sam_probs) == 1)
      }
      bin_out <- out %% N_Y_BIN
      bin_out[bin_out == 0] <- N_Y_BIN
      add_mat <- t(sapply(bin_out, function(x) (1:N_Y_BIN) == x) / control$N_ITS)
      prob_mat <- prob_mat + add_mat
      out_mat <- cbind(out_mat, bin_out)
    }
    
    y_sam_mat <- out_mat
    y_sam_mat[] <- runif(length(y_sam_mat), min = y_bin_df$l[out_mat], max = y_bin_df$u[out_mat])
    Sig_subj <- cov(t(y_sam_mat))
    mu_subj <- rowMeans(y_sam_mat)
    log_Sig_subj <- cov(t(log(y_sam_mat)))
    log_mu_subj <- rowMeans(log(y_sam_mat))
    
    list_out[[subj_test]] <- list(mu = mu_subj, Sig = Sig_subj, mu_log = log_mu_subj, Sig_log = log_Sig_subj)
    
    if (control$do_plots) {
      ##############################################################################
      # Plot posterior summaries for this individual
      ##############################################################################
      # pdf(file = paste0("plots/temp.pdf"), width = 9, height = 9)
      out_y_mat<-out_y_mat_l<-out_y_mat_u <- out_mat
      out_y_mat[] <- y_bin_df$m[out_mat]
      out_y_mat_l[] <- y_bin_df$l[out_mat]
      out_y_mat_u[] <- y_bin_df$u[out_mat]
      dat_all$age_yrs <- dat_all$age_days / 365.25
      y_qs <- cbind(t(apply(out_y_mat_l, 1, function(v) quantile(v, c(0.05, .25)))),
                    (apply(out_y_mat, 1, function(v) quantile(v, c(.5)))),
                    t(apply(out_y_mat_u, 1, function(v) quantile(v, c(0.75, 0.95)))))
      age_bin_df$m <- (age_bin_df$l + age_bin_df$u) / 2
      # matplot(x = age_bin_df$m, y_qs, 
      #         lty = c(3, 2, 1, 2, 3), ty = "l", col = 1, lwd = c(1, 1, 2, 1, 1),
      #         ylim = c(15, 45), main = subj_test)
      # points(dat_all$age_yrs, dat_all$value)
      plot(x = age_bin_df$m, y = mu_subj, ty = "l", ylim = c(15, 45), main = subj_test)
      lines(x = age_bin_df$m, y = mu_subj + qnorm(.975) * sqrt(diag(Sig_subj)), lty = 2)
      lines(x = age_bin_df$m, y = mu_subj - qnorm(.975) * sqrt(diag(Sig_subj)), lty = 2)
      lines(x = age_bin_df$m, y = mu_subj + qnorm(.75) * sqrt(diag(Sig_subj)), lty = 3)
      lines(x = age_bin_df$m, y = mu_subj - qnorm(.75) * sqrt(diag(Sig_subj)), lty = 3)
      lines(x = age_bin_df$m, y = y_mn, ty = "l", ylim = c(15, 45), main = subj_test, col = 2)
      lines(x = age_bin_df$m, y = y_mn + qnorm(.975) * sqrt(diag(y_var)), lty = 2, col = 2)
      lines(x = age_bin_df$m, y = y_mn - qnorm(.975) * sqrt(diag(y_var)), lty = 2, col = 2)
      n_sam_plot <- 100
      abline(v = AGE_MN_TRAIN + c(-1, 1) * AGE_RADIUS_TRAIN)
      sams <- MASS::mvrnorm(n = n_sam_plot, mu = mu_subj, 
                            Sigma = Sig_subj) 
      points(dat_all$age_yrs, dat_all$value)
      # for (j in 1:n_sam_plot) {
      #   lines(sams[j,  ], lty = 1, col = 4)
      # }
      # dev.off()
    }
  }
}

if (control$do_plots) {
  dev.off()
}


dir.create("output/mv_gauss_posteriors", showWarnings = FALSE)
saveRDS(list_out, file = paste0("output/mv_gauss_posteriors/mv_gauss_posterior_task_", TASK_ID, ".RDS"))

