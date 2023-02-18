
#' Function to calculate age bins
#' @param control List of control parameters
get_age_bin_df <- function(control) {
  age_bin_df <- data.frame(l = seq(control$MIN_AGE, control$MAX_AGE, by = control$AGE_BIN_WID))
  age_bin_df$u <- age_bin_df$l + control$AGE_BIN_WID
  age_bin_df$bin <- 1:nrow(age_bin_df)
  age_bin_df$bin_name <- paste0("[", age_bin_df$l, ", ", age_bin_df$u, ")")
  age_bin_df$low_day <- round(age_bin_df$l * 365.25)
  age_bin_df$upp_day <- round(age_bin_df$u * 365.25) - 1
  return(age_bin_df)
}
