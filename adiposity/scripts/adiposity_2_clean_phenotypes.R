# Author: Samvida S. Venkatesh
# Date: 13/02/21

library(tidyverse)
library(lubridate)

# Read and subset required data ----

# Adiposity phenotypes
adiposity <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/raw_adiposity.rds")
PHENOTYPES <- names(adiposity)

# Covariates file 
general_covars <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/QCd_covariates.txt",
                             sep = "\t", header = T, comment.char = "$",
                             stringsAsFactors = F)

# Bariatric surgery codes
bariatric_codes <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/bariatric_surgery_records.txt",
                              sep = "\t", header = T, comment.char = "$",
                              stringsAsFactors = F)
bariatric_codes$event_dt <- as.Date(bariatric_codes$event_dt, "%Y-%m-%d")

# # Pregnancy codes
# pregnancy_codes <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/pregnancy_records.txt",
#                               sep = "\t", header = T, comment.char = "$",
#                               stringsAsFactors = F)
# pregnancy_codes$start <- as.Date(pregnancy_codes$start, "%Y-%m-%d")
# pregnancy_codes$end <- as.Date(pregnancy_codes$end, "%Y-%m-%d")
# pregnancy_codes$preg_interval <- interval(pregnancy_codes$start, 
#                                           pregnancy_codes$end)
# pregnancy_intervals <- split(pregnancy_codes[, "preg_interval"],
#                              pregnancy_codes$eid)

# Convert weight measures to BMI and vice-versa ----

# weight to BMI
weight_to_BMI <- adiposity$weight
weight_to_BMI$height <- general_covars[match(weight_to_BMI$eid, 
                                             general_covars$eid), "height"]
weight_to_BMI <- weight_to_BMI %>% mutate(weight = value,
                                          BMI = value / (height/100)^2)
weight_to_BMI$value <- weight_to_BMI$BMI
weight_to_BMI <- weight_to_BMI[complete.cases(weight_to_BMI$value), ]

weight_to_BMI <- bind_rows(adiposity$BMI, 
                           weight_to_BMI[, colnames(adiposity$BMI)]) %>%
  arrange(eid, event_dt)
adiposity$BMI <- weight_to_BMI

# BMI to weight
BMI_to_weight <- adiposity$BMI
BMI_to_weight$height <- general_covars[match(BMI_to_weight$eid, 
                                             general_covars$eid), "height"]
BMI_to_weight <- BMI_to_weight %>% mutate(BMI = value,
                                          weight = value * (height/100)^2)
BMI_to_weight$value <- BMI_to_weight$weight
BMI_to_weight <- BMI_to_weight[complete.cases(BMI_to_weight$value), ]

BMI_to_weight <- bind_rows(adiposity$weight, 
                           BMI_to_weight[, colnames(adiposity$weight)]) %>%
  arrange(eid, event_dt)
adiposity$weight <- BMI_to_weight

# Functions to clean data ----

## Age ----

ageFilter <- function (df, qc_log_file) {
  cleaned <- subset(df, df$age_event >= 20)
  # Report QC metrics
  sink(qc_log_file, append = T)
  cat(paste0("**FILTER** EXCLUDED, Age at measurement < 20 years: ", "\n",
             "\t", "Number of measurements = ", 
             dim(df)[1] - dim(cleaned)[1], "\n",
             "\t", "Number of individuals = ", 
             length(unique(df$eid)) - length(unique(cleaned$eid)), "\n"))
  sink()
  return (cleaned)
}

## Bariatric surgery ----

barSurgeryFilter <- function (df, qc_log_file) {
  
  # Exclude individuals with history of bariatric surgery
  bar_hist_EIDS <- bariatric_codes$eid[is.na(bariatric_codes$event_dt)]
  # Record date of surgery for individuals with that information
  df$post_bariatric_surgery_flag <- 
    bariatric_codes$event_dt[match(df$eid, bariatric_codes$eid)]
  df$post_bariatric_surgery_flag <- df$event_dt > 
    df$post_bariatric_surgery_flag | df$eid %in% bar_hist_EIDS
  df$post_bariatric_surgery_flag[is.na(df$post_bariatric_surgery_flag)] <- F
  
  cleaned <- subset(df, !df$post_bariatric_surgery_flag)
  
  # Report QC metrics
  sink(qc_log_file, append = T)
  cat(paste0("**FILTER** EXCLUDED, Measured post-bariatric surgery: ", "\n",
             "\t", "Number of measurements = ", 
             sum(df$post_bariatric_surgery_flag, na.rm = T), "\n",
             "\t", "Number of individuals = ", 
             length(unique(df$eid)) - length(unique(cleaned$eid)), "\n"))
  sink()
  
  return (cleaned)
  
}

# ## Pregnancy ----
# 
# pregFilter <- function (df, qc_log_file) {
#   
#   no_preg <- subset(df, !df$eid %in% pregnancy_codes$eid)
#   possible_preg <- subset(df, df$eid %in% pregnancy_codes$eid)
#   
#   # Record pregnancy interval
#   possible_preg$eid <- as.character(possible_preg$eid)
#   # Adapted code from: https://github.com/tidyverse/lubridate/issues/658
#   possible_preg$preg_flag <- 
#     unlist(lapply(sapply(1:dim(possible_preg)[1], 
#                          function(i) possible_preg$event_dt[i] %within% 
#                            pregnancy_intervals[possible_preg$eid[i]]), 
#                   function (x) any(x)))
#   
#   # Combine 
#   no_preg$preg_flag <- F
#   cleaned <- bind_rows(no_preg, possible_preg) %>% arrange(eid, event_dt)
#   
#   indivs_flagged <- unique(cleaned$eid[which(cleaned$preg_flag)])
#   
#   # Report QC metrics
#   sink(qc_log_file, append = T)
#   cat(paste0("**FILTER** FLAGGED, Measured during pregnancy: ", "\n",
#              "\t", "Number of measurements = ", 
#              sum(cleaned$preg_flag, na.rm = T), "\n",
#              "\t", "Number of individuals = ", 
#              length(indivs_flagged), "\n"))
#   sink()
#   
#   return (cleaned)
#   
# }

## Implausible values ----

IMPLAUSIBLE_MINS <- c(15, 40, 10, 0.5)
IMPLAUSIBLE_MAXES <- c(200, 160, 70, 1.4)
names(IMPLAUSIBLE_MINS) <- PHENOTYPES
names(IMPLAUSIBLE_MAXES) <- PHENOTYPES

implausibleFilter <- function (df, p, qc_log_file) {
  # Remove implausible values defined above
  remove_mins <- df[df$value <= IMPLAUSIBLE_MINS[p], ]
  remove_maxes <- df[df$value >= IMPLAUSIBLE_MAXES[p], ]
  
  cleaned <- subset(df, df$value > IMPLAUSIBLE_MINS[p] & 
                      df$value < IMPLAUSIBLE_MAXES[p])
  
  # Report QC metrics
  sink(qc_log_file, append = T)
  cat(paste0("**FILTER** EXCLUDED, Implausible value below ", 
             IMPLAUSIBLE_MINS[p], ":", "\n",
             "\t", "Number of measurements = ", 
             dim(remove_mins)[1], "\n",
             "\t", "Number of individuals = ", 
             sum(!unique(remove_mins$eid) %in% unique(cleaned$eid)), "\n",
             "**FILTER** EXCLUDED, Implausible value above ", 
             IMPLAUSIBLE_MAXES[p], ":", "\n",
             "\t", "Number of measurements = ", 
             dim(remove_maxes)[1], "\n",
             "\t", "Number of individuals = ", 
             sum(!unique(remove_maxes$eid) %in% unique(cleaned$eid)), "\n"))
  sink()
  return (cleaned)
}

## Extreme values ----

extremeFilter <- function (df, qc_log_file) {
  # Remove values +/- 5 S.D. away from the population mean
  popn_mean <- mean(df$value, na.rm = T)
  popn_sd <- sd(df$value, na.rm = T)
  max_outlier <- popn_mean + 5*popn_sd
  min_outlier <- popn_mean - 5*popn_sd
  cleaned <- subset(df, df$value <= max_outlier & df$value >= min_outlier)
  # Report QC metrics
  sink(qc_log_file, append = T)
  cat(paste0("**FILTER** EXCLUDED, Value > 5 S.D. away from population mean: ", 
             "\n",
             "\t", "Number of measurements = ", 
             dim(df)[1] - dim(cleaned)[1], "\n",
             "\t", "Number of individuals = ", 
             length(unique(df$eid)) - length(unique(cleaned$eid)), "\n"))
  sink()
  return (cleaned)
}

## Timepoints (longitudinal) ----

longitFilter <- function (df, qc_log_file) {
  # Remove individuals with fewer than 2 measurements post-cleaning
  nmeasures <- df %>% group_by(eid) %>% summarise(n = n())
  remove_ids <- nmeasures$eid[nmeasures$n < 2]
  cleaned <- subset(df, !df$eid %in% remove_ids)
  # Report QC metrics
  sink(qc_log_file, append = T)
  cat(paste0("**FILTER** EXCLUDED, Only has 1 post-QC measurement: ", 
             "\n",
             "\t", "Number of individuals = ", 
             length(remove_ids), "\n"))
  sink()
  return (cleaned)
}

# Apply QC filters ----

cleaned_adiposity <- lapply(PHENOTYPES, function (p) {
  
  log_file_p <- paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/log_files/qc_log_file_",
                       p, ".txt")
  cleaned <- ageFilter(adiposity[[p]], log_file_p)
  cleaned <- barSurgeryFilter(cleaned, log_file_p)
  #cleaned <- pregFilter(cleaned, log_file_p)
  cleaned <- implausibleFilter(cleaned, p, log_file_p)
  cleaned <- extremeFilter(cleaned, log_file_p)
  cleaned <- longitFilter(cleaned, log_file_p)
  
  return (cleaned)
  
})
names(cleaned_adiposity) <- PHENOTYPES

saveRDS(cleaned_adiposity, "/well/lindgren/UKBIOBANK/samvida/adiposity/QCd_adiposity.rds")