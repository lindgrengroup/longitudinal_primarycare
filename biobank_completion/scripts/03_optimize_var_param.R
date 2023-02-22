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
  TASK_ID <- 39
  N_TASKS <- 100
}

all_traj <- readRDS(file = "output/all_traj.RDS")
glob_mn <- Reduce(f = '+', x = lapply(all_traj, function(x) x$mu_log)) / length(all_traj)

set.seed(TASK_ID)
subj_opt <- sample(names(all_traj), size = control$N_SUBJ_PER_OPTIM)
inds_use <- 1:n_age
obj <- function(par) {
  Sig_eps <- Sig_eps_fn(len = exp(par[1]), 
                         sig_GP = exp(par[2]), 
                         sig_RI = exp(par[3]), 
                         sig_resid = exp(par[4]))
  densv <<- c()
  log_dens <- 0
  for (j in subj_opt) {
    curr_dens <- mvtnorm::dmvnorm(x = all_traj[[j]]$mu_log[inds_use], mean = glob_mn[inds_use], 
                     sigma = (all_traj[[j]]$Sig_log + Sig_eps)[inds_use, inds_use], log = TRUE)
    densv[j] <<- curr_dens
    log_dens <- log_dens + curr_dens
  }
  print(c(len = exp(par[1]),
              sig_GP = exp(par[2]),
              sig_RI = exp(par[3]),
              sig_resid = exp(par[4])))
  return(-log_dens)
}

par_start <- log(c(1, 1, 1, 1))
out <- optim(par = par_start, fn = obj, control = list(trace = 10))
dir.create("output/par_optim", showWarnings = FALSE)
saveRDS(out, file = paste0("output/par_optim/optim_out_", TASK_ID, ".RDS"))

# 
# mu_plot <- NULL
# for (j in subj_opt) {
#   mu_plot <- cbind(mu_plot, all_traj[[j]]$mu_log[inds_use])
# }
# sig_plot <- NULL
# for (j in subj_opt) {
#   sig_plot <- cbind(sig_plot, diag(all_traj[[j]]$Sig_log[inds_use, inds_use]))
# }
# 
# pdf("plots/subj_outlier.pdf", 12, 8)
# matplot(mu_plot, type = "l")
# dev.off()
# 
# pdf("plots/sig_subj_outlier.pdf", 12, 8)
# matplot(sig_plot, type = "l")
# dev.off()
