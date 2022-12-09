# Author: Samvida S. Venkatesh
# Date: 04/10/21

library(tidyverse)
library(lubridate)

# Read data ----

# Cleaned phenotypes from full primary care and UKB data 
dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/gp_main_data_passed_longit_filter.rds")
PHENOTYPES <- names(dat)

# Start QC log file
qc_log <- "/well/lindgren/UKBIOBANK/samvida/full_primary_care/qc/corrected_indiv_qc.txt"

# Filter for longitudinal data ----

longitFilter <- function (df, p) {
  
  # Remove multiple measurements at same timepoint 
  # (remove the one farthest from median)
  cleaned <- df %>% group_by(eid) %>% 
    mutate(median_value = median(value, na.rm = T),
           dist_to_med = abs(value - median_value)) 
  cleaned <- cleaned %>% group_by(eid, event_dt) %>% 
    mutate(remove = n() > 1 & dist_to_med == max(dist_to_med))
  cleaned <- subset(cleaned, !cleaned$remove)
  
  # Remaining values should be identical, so simply take 
  # the first
  cleaned <- cleaned %>% distinct(eid, data_provider, event_dt, 
                                  age_event, .keep_all = T)
  
  # Remove individuals with fewer than 2 measurements post-cleaning
  nmeasures <- cleaned %>% group_by(eid) %>% summarise(n = n())
  remove_ids <- nmeasures$eid[nmeasures$n < 2]
  cleaned <- subset(cleaned, !cleaned$eid %in% remove_ids)
  # Report QC metrics
  sink(qc_log, append = T)
  cat(paste0("** PHENOTYPE **", p, "\n",  
             "**FILTER** EXCLUDED, Only has 1 post-QC measurement: ", 
             "\n",
             "\t", "Number of individuals = ", 
             length(remove_ids), "\n"))
  sink()
  
  return (cleaned)
}

# Apply QC filters ----

cleaned_dat <- lapply(PHENOTYPES, function (p) {
  cleaned <- longitFilter(dat[[p]], p)
  return (cleaned)
})
names(cleaned_dat) <- PHENOTYPES

# Remove data-points causing large jumps ----

# Calculate jump as time-adjusted fold-change
jumps <- lapply(cleaned_dat, function (df) {
  
  # If two measurements are taken at the same age, 
  # remove the one farther from the overall median value
  dups <- df %>% group_by(eid) %>%
    mutate(dist_to_med = abs(value - median(value)),
           # check if there are two values at the same age
           dup_marker = duplicated(age_event) | duplicated(age_event, fromLast = T))
  
  res <- subset(dups, !dups$dup_marker)
  if (any(dups$dup_marker)) {
    no_dups <- subset(dups, dups$dup_marker) %>% group_by(eid, age_event) %>%
      mutate(remove = dist_to_med == max(dist_to_med))
    no_dups <- subset(no_dups, !no_dups$remove)
    res <- bind_rows(res, no_dups)
  }
  
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
})

indiv_cleaned <- lapply(jumps, function (df) {
  
  if (any(df$extreme, na.rm = T)) {
    extremes <- df[df$eid %in% df$eid[df$extreme], ]
    # Remove the point farthest from the median (that causes the jump)
    extremes <- extremes %>% group_by(eid) %>%
      mutate(median_value = median(value),
             dist_to_med = abs(value - median_value),
             remove = dist_to_med == max(dist_to_med))
    extremes <- subset(extremes, !extremes$remove)
    cleaned <- bind_rows(df[!df$eid %in% df$eid[df$extreme], ],
                         extremes[, colnames(df)])
  } else cleaned <- df
  
  # Save relevant information
  cleaned <- cleaned[, c("eid", "data_provider", "event_dt",
                         "age_event", "value")]
  # Ensure eids are saved as characters (not factors)
  cleaned$eid <- as.character(cleaned$eid)
  
  return (cleaned)
})

# Log number of observations and individuals removed by jumps ----

lapply(PHENOTYPES, function (p) {
  orig_obs <- dim(cleaned_dat[[p]])[1]
  orig_indivs <- length(unique(cleaned_dat[[p]]$eid))
  
  cleaned_obs <- dim(indiv_cleaned[[p]])[1]
  cleaned_indivs <- length(unique(indiv_cleaned[[p]]$eid))
  
  sink(qc_log, append = T)
  cat(paste0("** PHENOTYPE **", p, "\n",  
             "**FILTER** EXCLUDED, Unreasonable jump > 5 S.D. away from popn.: ", 
             "\n",
             "\t", "Number of measurements = ", 
             orig_obs - cleaned_obs, "\n",
             "\t", "Number of individuals = ", 
             orig_indivs - cleaned_indivs, "\n"))
  sink()
})

saveRDS(indiv_cleaned, 
        "/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/corrected_indiv_qcd_data.rds")