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

all_traj <- readRDS(file = "output/all_traj.RDS")
par_mat <- NULL
for (TASK_ID in 1:control$NTASKS_03) {
  out <- readRDS(file = paste0("output/par_optim/optim_out_", TASK_ID, ".RDS"))
  par_mat <- rbind(par_mat, t(out$par))
}

par_use <- colMeans(exp(par_mat))
Sig_eps <- Sig_eps_fn(len = par_use[1], sig_GP = par_use[2], sig_RI = par_use[3], sig_resid = par_use[4])
subj_split <- subj_split_fn(all_subj, N_TASKS = N_TASKS)  
subj_do <- subj_split[[TASK_ID]]

vec_combine<-mat_combine <- 0
for (subj_curr in subj_do) {
  print(subj_curr)
  meta_curr <- meta[meta$subj_id == subj_curr, ]
  x_curr <- c(meta_curr$sex, meta_curr[snp_names] == 1, meta_curr[snp_names] == 2)
  prec_curr <- solve(Sig_eps + all_traj[[j]]$Sig)
  mat_combine <- mat_combine + kronecker(X = x_curr %*% t(x_curr), Y = prec_curr)
  vec_combine <- vec_combine + kronecker(X = x_curr, Y = prec_curr %*% all_traj[[j]]$mu)
  gc()
}

dir.create("output/linear_model_fitting", showWarnings = FALSE)
saveRDS(list(mat = mat_combine, vec = vec_combine), file = paste0("output/linear_model_fitting/mat_vec_combine_task_", TASK_ID, ".RDS"))


