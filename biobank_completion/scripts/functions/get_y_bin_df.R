
#' Function to calculate trait bins
#' @param control List of control parameters
#' @param y_vec Vector of all trait values (over which to calculate quantile bins) 
get_y_bin_df <- function(control, y_vec) {
  bin_lower <- quantile(y_vec, probs = seq(0, 1 - 1 / control$N_Y_BIN, by = 1 / control$N_Y_BIN))
  bin_upper <- quantile(y_vec, probs = seq(1 / control$N_Y_BIN, 1, by = 1 / control$N_Y_BIN))
  y_bin_df <- data.frame(l = bin_lower, u = bin_upper)
  y_bin_df$m <- (y_bin_df$l + y_bin_df$u) / 2
  return(y_bin_df)
}
