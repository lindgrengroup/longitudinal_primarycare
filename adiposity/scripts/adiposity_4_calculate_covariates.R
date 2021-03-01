# Author: Samvida S. Venkatesh
# Date: 19/02/21

library(tidyverse)
library(lubridate)

# Read files ----

adiposity <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/QCd_adiposity.rds")
general_covars <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/QCd_covariates.txt",
                             sep = "\t", header = T, comment.char = "$",
                             stringsAsFactors = F)
PHENOTYPES <- names(adiposity)

# Design covariate files for each adiposity trait ----

NPCs = 21

# Get BMI data to calculate baseline
QCd_for_baseline <- adiposity$BMI
QCd_for_baseline <- QCd_for_baseline %>% arrange(eid, event_dt)
QCd_for_baseline <- split(QCd_for_baseline[, c("event_dt", "type", "value")], 
                          f = QCd_for_baseline$eid)

QCd_covars <- lapply(PHENOTYPES, function (p) { 
  df <- adiposity[[p]]
  # Calculate baseline and follow-up years covariates
  calc_covars <- df %>% group_by(eid) %>% 
    arrange(eid, event_dt) %>%
    summarise(baseline_date = first(event_dt),
              baseline_age = first(age_event),
              age_sq = baseline_age^2,
              FUyrs = interval(first(event_dt), last(event_dt)) / years(1))
  
  # Get baseline BMI from file (nearest BMI measure to first adiposity)
  calc_covars$baseline_BMI <- sapply(1:dim(calc_covars)[1], function (i) {
    dat_sub <- QCd_for_baseline[[calc_covars$eid[i]]]
    if (!is.null(dat_sub)) {
      closest_measure <- which.min(abs(calc_covars$baseline_date[i] - 
                                         dat_sub$event_dt))
      res <- dat_sub$value[closest_measure]
    } else res <- NA
    return (res)
  })
  
  # Merge with previously calculated covariates
  cleaned <- merge(calc_covars, general_covars, by = "eid")
  
  # Remove observations missing baseline BMI
  sink(paste0("log_files/covariate_QC_", p, ".txt"), append = T)
  cat(paste0("**FILTER** EXCLUDED, No baseline BMI: ", 
             length(which(is.na(cleaned$baseline_BMI))), "\n"))
  sink()
  cleaned <- subset(cleaned, !is.na(cleaned$baseline_BMI))
  
  # Set missing and inconsistent ancestry to "other"
  cleaned$ancestry <- ifelse(cleaned$ancestry == "missing or inconsistent",
                             "other", cleaned$ancestry)
  
  sink(paste0("log_files/covariate_QC_", p, ".txt"), append = T)
  cat(paste0("**FILTER** EXCLUDED, General QC failed
             (genotyping, relatedness, recommended exclusions): ", "\n",
             dim(calc_covars)[1] - dim(cleaned)[1], "\n"))
  sink()
  
  # Remove individuals without a height measurement
  sink(paste0("log_files/covariate_QC_", p, ".txt"), append = T)
  cat(paste0("**FILTER** EXCLUDED, No height: ", 
             length(which(is.na(cleaned$height))), "\n"))
  sink()
  cleaned <- subset(cleaned, !is.na(cleaned$height))
  
  # Remove individuals missing any other covariate
  res <- cleaned[, c("eid", "sex", "ancestry",
                     "baseline_age", "age_sq", 
                     "height", "baseline_BMI", "FUyrs", 
                     "genotyping_array", paste0("PC", 1:NPCs))]
  res <- res[complete.cases(res), ]
  sink(paste0("log_files/covariate_QC_", p, ".txt"), append = T)
  cat(paste0("**FILTER** EXCLUDED, Missing any other covariate: ", 
             dim(cleaned)[1] - dim(res)[1], "\n"))
  sink()
  
  return(res)
})

names(QCd_covars) <- PHENOTYPES

# Stratify on sex and ancestry ----

for (p in PHENOTYPES) {
  sink(paste0("log_files/stratified_counts_", p, ".txt"), append = T)
  cat(paste0("Number of individuals: ", "\n"))
  print(with(QCd_covars[[p]], table(ancestry, sex)))
  cat("\n")
  sink()
}

# Save
saveRDS(QCd_covars, "/well/lindgren/UKBIOBANK/samvida/adiposity/adiposity_covars.rds")