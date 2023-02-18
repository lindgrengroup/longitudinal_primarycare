source("scripts/00_load_data.R")
##############################################################################
# Combine MV Gaussian trajectories
##############################################################################
all_traj <- list()
for (TASK_ID in 1:control$NTASKS_02) {
  print(TASK_ID)
  task_traj <- readRDS(file = paste0("output/mv_gauss_posteriors/mv_gauss_posterior_task_", TASK_ID, ".RDS"))
  all_traj <- c(all_traj, task_traj)
}
saveRDS(all_traj, file = "output/all_traj.RDS")
