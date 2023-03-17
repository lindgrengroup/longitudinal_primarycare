#' Function to split subjects into subsets
#' @param all_subj Vector of all subjects
subj_split_fn <- function(all_subj, N_TASKS) {
  n_subj <- length(all_subj)
  subj_seq_ind <- floor(seq(1, n_subj, length.out = N_TASKS + 1))
  subj_split <- sapply(1:N_TASKS, function(j) all_subj[subj_seq_ind[j]:subj_seq_ind[j + 1]])
  return(subj_split)
}
