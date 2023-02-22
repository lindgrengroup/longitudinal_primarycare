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
  N_TASKS <- 50
}

mat_combine <- vec_combine <- 0
n_to_combine <- control$NTASKS_04 / control$NTASKS_04a
tasks_to_combine <- (TASK_ID - 1) * n_to_combine + 1:n_to_combine
for (TASK04_ID in tasks_to_combine) {
  res_list_curr <- readRDS(file = paste0("output/linear_model_fitting/mat_vec_combine_task_", TASK04_ID, ".RDS"))
  mat_combine <- mat_combine + res_list_curr$mat
  vec_combine <- vec_combine + res_list_curr$vec
}
saveRDS(list(mat = mat_combine, vec = vec_combine), file = paste0("output/linear_model_fitting/mat_vec_combine_task04a_", TASK_ID, ".RDS"))
