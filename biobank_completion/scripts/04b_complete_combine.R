source("scripts/00_load_data.R")

mat_combine <- vec_combine <- 0
for (TASK04a_ID in 1:control$NTASKS_04a) {
  res_list_curr <- readRDS(file = paste0("output/linear_model_fitting/mat_vec_combine_task04a_", TASK04a_ID, ".RDS"))
  mat_combine <- mat_combine + res_list_curr$mat
  vec_combine <- vec_combine + res_list_curr$vec
}
saveRDS(list(mat = mat_combine, vec = vec_combine), file = paste0("output/linear_model_fitting/mat_vec_combine_task04a_", TASK_ID, ".RDS"))