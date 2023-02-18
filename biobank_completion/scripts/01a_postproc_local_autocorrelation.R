source("scripts/00_load_data.R")
TASK_ID <- 1
N_TASKS <- 1000
joint_dist_mat <- matrix(0, N_Y_BIN, N_Y_BIN)
joint_dist_mat_multi_list <- list()
for (j in 1:control$N_GRAD_CLASSES) {
  joint_dist_mat_multi_list[[j]] <- rep(list(joint_dist_mat), N_AGE_BIN)
}
for (TASK_ID in 1:N_TASKS) {
  joint_dist_mat_multi_list_curr <- readRDS(file = paste0("output/joint_dist_mat/joint_dist_mat_multi_list_task_", TASK_ID, ".RDS"))
  for (j in 1:control$N_GRAD_CLASSES) {
    for (age_bin_curr in 1:N_AGE_BIN) {
      joint_dist_mat_multi_list[[j]][[age_bin_curr]] <- joint_dist_mat_multi_list[[j]][[age_bin_curr]] + joint_dist_mat_multi_list_curr[[j]][[age_bin_curr]]
    }
  }
}

##############################################################################
# Post process autocorrelation structure 
##############################################################################
joint_dist_mat_multi_list_prob_dist <- lapply(joint_dist_mat_multi_list, function(curr_list) {
  lapply(curr_list, function(M){
    M <- M / sum(M)
  })
})

smooth_trans_mats <- TRUE
if (smooth_trans_mats) {
  joint_dist_mat_multi_list_prob_dist_smoothed <- lapply(joint_dist_mat_multi_list_prob_dist, function(curr_list) {
    curr_list_out <- curr_list
    for (i in 3:(length(curr_list) - 2)) {
      curr_list_out[[i]] <- (curr_list[[i - 2]] + curr_list[[i - 1]] + curr_list[[i]] + curr_list[[i + 1]] + curr_list[[i + 2]]) / 5
      curr_list_out[[i]][is.na(curr_list_out[[i]])] <- curr_list[[i]][is.na(curr_list_out[[i]])]
    }
    return(curr_list_out)
  })
  
  joint_dist_mat_multi_list_conditional <- lapply(joint_dist_mat_multi_list_prob_dist_smoothed, function(curr_list) {
    lapply(curr_list, function(M){
      M <- M / rowSums(M)
    })
  })
} else {
  joint_dist_mat_multi_list_conditional <- lapply(joint_dist_mat_multi_list_prob_dist, function(curr_list) {
    lapply(curr_list, function(M){
      M <- M / rowSums(M)
    })
  })
}

##############################################################################
# Plot transition matrices
##############################################################################
if (control$do_plots) {
  pdf(file = "plots/transmats.pdf", width = 9, height = 12)
  par(mfrow = c(4, control$N_GRAD_CLASSES), mar = c(1, 1, 1, 1))
  for (i in 1:(N_AGE_BIN - 1)) {
    for (j in 1:control$N_GRAD_CLASSES) {
      image(joint_dist_mat_multi_list_conditional[[j]][[i]])
      abline(0, 1)
    }
  }
  dev.off()
}

##############################################################################
# Plot differenced transition matrices
##############################################################################
if (control$do_plots) {
  pdf(file = "plots/transmats_diff.pdf", width = 9, height = 12)
  par(mfrow = c(10, control$N_GRAD_CLASSES), mar = c(2, 1, 2, 1), oma = c(2, 4, 2, 4))
  cexax <- .5
  for (i in 1:(N_AGE_BIN - 1)) {
    for (j in 1:control$N_GRAD_CLASSES) {
      M <- joint_dist_mat_multi_list_conditional[[j]][[i]]
      max_n_diff <- 10
      M_diff <- matrix(0, N_Y_BIN, 2 * max_n_diff + 1)
      for (from_ind in 1:N_Y_BIN) {
        for (diff_curr in (-max_n_diff):max_n_diff) {
          if (from_ind + diff_curr <= N_Y_BIN & from_ind + diff_curr >= 1) {
            M_diff[from_ind, diff_curr + max_n_diff + 1] <- M[from_ind, from_ind + diff_curr]
          }
        }
      }
      image(x = 1:N_Y_BIN, y = (-max_n_diff):max_n_diff, z = M_diff, ylab = "", cex.axis = .7, las = 1)
      abline(h = 0)
      legend(x = "topleft", legend = paste0("Previous change = ", j - 3), cex = .6)
      trans_text <- paste0("From ", age_bin_df[i, "bin_name"], " to ", age_bin_df[i + 1, "bin_name"])
      if (j == 1) {
        mtext(side = 2, text = "Change in y bin #", line = 2, cex = cexax)
        mtext(side = 2, text = trans_text, line = 3, cex = cexax)
      }
      mtext(side = 1, text = "Originating y bin #", line = 2, cex = cexax)
    }
  }
  dev.off()
}

##############################################################################
# Create large transition matrix including location and gradient
##############################################################################
raw_trans <- joint_dist_mat_multi_list_conditional
trans_mats <- list()
for (age_bin_num in 1:(N_AGE_BIN - 1)) {
  trans_mat_curr <- NULL
  for (grad_from in 1:control$N_GRAD_CLASSES) {
    add_wide_mat <- NULL
    for (grad_to in 1:control$N_GRAD_CLASSES) {
      add_mat <- raw_trans[[grad_from]][[age_bin_num]]
      # Repairing rows lacking data with identity transform
      rows_to_repair <- which(rowMeans(is.na(add_mat)) > 0)
      add_mat[rows_to_repair, ] <- 0
      add_mat[cbind(rows_to_repair, rows_to_repair)] <- 1
      max_abs_grad <- (control$N_GRAD_CLASSES - 1) / 2
      if (grad_to == 1) {
        add_mat[col(add_mat) - row(add_mat) > -max_abs_grad] <- 0
      } else if (grad_to == control$N_GRAD_CLASSES) {
        add_mat[col(add_mat) - row(add_mat) < max_abs_grad] <- 0
      } else {
        mid_class <- mean(1:control$N_GRAD_CLASSES)
        add_mat[col(add_mat) - row(add_mat) != grad_to - mid_class] <- 0
      }
      add_wide_mat <- cbind(add_wide_mat, add_mat)
    }
    trans_mat_curr <- rbind(trans_mat_curr, add_wide_mat)
  }
  trans_mats[[age_bin_num]] <- trans_mat_curr
}
saveRDS(trans_mats, file = "output/trans_mats.RDS")

##############################################################################
# Check effect of gradient
##############################################################################
if (control$do_plots) {
  for (from_class in 1:control$N_GRAD_CLASSES) {
    print(from_class)
    for (age_bin in 1:N_AGE_BIN) {
      print(c(sum(raw_trans[[from_class]][[age_bin]][lower.tri(raw_trans[[from_class]][[age_bin]])], na.rm = T),
              sum(diag(raw_trans[[from_class]][[age_bin]]), na.rm = T),
              sum(raw_trans[[from_class]][[age_bin]][upper.tri(raw_trans[[from_class]][[age_bin]])], na.rm = T)))
    }
  }
}

##############################################################################
# Plot some full transition matrices
##############################################################################
if (control$do_plots) {
  pdf(file = "plots/full_trans.pdf", width = 9, height = 9)
  for (age_bin_num in 1:10) {
    image(trans_mats[[age_bin_num]])
  }
  dev.off()
}




