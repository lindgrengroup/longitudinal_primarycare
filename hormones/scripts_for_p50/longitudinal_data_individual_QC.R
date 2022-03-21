# Author: Samvida S. Venkatesh
# Date: 04/10/21

library(tidyverse)
library(lubridate)

# Read data ----

# Data that has passed population-level QC
# Must have fields for: 
# 1. individual id (eid), 
# 2. hormone value (value), and 
# 3. age at hormone measurement (age_event) 
# Optional:
# 4. date of hormone measurement (event_dt)

dat <- readRDS("data_passed_longit_filter.rds")

# Start QC log file
qc_log <- "indiv_qc_logs.txt"

# General purpose function to calculate distance from individual's 
# median value for any given dataframe

getDistToMedian <- function (df) {
  res <- df %>% group_by(eid) %>% 
    mutate(median_value = median(value, na.rm = T),
           dist_to_med = abs(value - median_value)) %>%
    ungroup()
  return (res)
}

# Filter for longitudinal data ----

# Remove multiple measurements at same timepoint 
# (remove the one farthest from median)
cleaned <- getDistToMedian(dat) %>%
  # Group by age_event if there is no event_date column
  group_by(eid, event_dt) %>% 
  mutate(remove = n() > 1 & dist_to_med == max(dist_to_med))

cleaned <- subset(cleaned, !cleaned$remove)

# Remaining values should be identical, so simply take 
# the first value
# Retain any additional fields of interest (ex. data provider, age, etc.)
cleaned <- cleaned %>% distinct(eid, data_provider, 
                                event_dt, age_event, 
                                .keep_all = T)

# Report QC metrics
sink(qc_log, append = T)
cat(paste0("**FILTER** EXCLUDED, Multiple measurements for same individual at same time-point: ", 
           "\n",
           "\t", "Number of measurements = ", 
           length(dat) - length(cleaned), "\n"))
sink()

# Remove data-points causing large jumps ----

calculateJump <- function (df) {
  # If two measurements are taken at the same age, 
  # remove the one farther from the overall median value
  # This should have been taken care of at the earlier step but this is a
  # sanity check
  dups <- getDistToMedian(df) %>%
    group_by(eid) %>%
    # check if there are two values at the same age
    mutate(dup_marker = duplicated(age_event) | duplicated(age_event, fromLast = T))
  
  res <- subset(dups, !dups$dup_marker)
  if (any(dups$dup_marker)) {
    no_dups <- subset(dups, dups$dup_marker) %>% 
      group_by(eid, age_event) %>%
      mutate(remove = dist_to_med == max(dist_to_med))
    no_dups <- subset(no_dups, !no_dups$remove)
    res <- bind_rows(res, no_dups)
  }
  
  # Calculate jump as time-adjusted fold-change
  res <- res %>% group_by(eid) %>% 
    # Ensure ordered
    arrange(age_event, .by_group = T) %>% 
    mutate(FC = abs(value - lag(value)) / lag(value),
           age_change = age_event - lag(age_event),
           jump = ifelse(FC == 0, NA, log2(FC/age_change)))
  
  # Mark jump as extreme if > 3 S.D. larger than mean jump
  mean_jump <- mean(res$jump, na.rm = T) 
  sd_jump <- sd(res$jump, na.rm = T)
  res$extreme <- res$jump > (mean_jump + (3*sd_jump))
  
  return (res)
}

# Add info on jump
jump_cleaned <- calculateJump(cleaned)

# Remove the point farthest from the median (that causes the jump)
if (any(jump_cleaned$extreme, na.rm = T)) {
  extremes <- subset(jump_cleaned, jump_cleaned$extreme)
  
  extremes <- getDistToMedian(extremes) %>%
    group_by(eid) %>%
    mutate(remove = dist_to_med == max(dist_to_med))
  extremes <- subset(extremes, !extremes$remove)
  
  # Add back the info on non-extreme points
  jump_cleaned <- bind_rows(subset(jump_cleaned, !jump_cleaned$extreme),
                            extremes[, colnames(jump_cleaned)])
} else res <- jump_cleaned

# Log number of observations and individuals removed by jumps ----

orig_obs <- dim(cleaned)[1]
orig_indivs <- length(unique(cleaned$eid))

cleaned_obs <- dim(jump_cleaned)[1]
cleaned_indivs <- length(unique(jump_cleaned$eid))

sink(qc_log, append = T)
cat(paste0("**FILTER** EXCLUDED, Unreasonable jump > 3 S.D. away from popn.: ", 
           "\n",
           "\t", "Number of measurements = ", 
           orig_obs - cleaned_obs, "\n",
           "\t", "Number of individuals = ", 
           orig_indivs - cleaned_indivs, "\n"))
sink()

# Save results ----

# Save relevant information
jump_cleaned <- jump_cleaned[, c("eid", "data_provider", "event_dt",
                                 "age_event", "value")]
# Ensure eids are saved as characters (not factors)
jump_cleaned$eid <- as.character(jump_cleaned$eid)

saveRDS(jump_cleaned, "indiv_qcd_data.rds")