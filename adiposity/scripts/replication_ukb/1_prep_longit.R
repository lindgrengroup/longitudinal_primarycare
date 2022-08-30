# Author: Samvida S. Venkatesh
# Date: 31/08/21

library(tidyverse)
library(lubridate)

# Read data ----

# UK Biobank main phenotype file
pheno <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "eid"

PHENOTYPES <- c("BMI", "Weight", "WC", "WHR", "Weight_change_1yr")
PHENO_CODES <- c("f.21001.*", "f.21002.*", "f.48.*", NA, "f.2306.*")
names(PHENO_CODES) <- PHENOTYPES

# Get relevant UKB data fo each phenotype ----

# Function to convert main pheno wide files to long format
wide_to_long <- function (wide_df, p) {
  res <- wide_df
  # Merge in centre and BMI codes to adj for BMI 
  extra_codes <- pheno[, c("eid", 
                           "f.54.0.0", "f.54.1.0", "f.54.2.0",
                           "f.21001.0.0", "f.21001.1.0", "f.21001.2.0")]
  res <- merge(res, extra_codes, by = "eid")
  # Pivot to long format
  # Special case for WHR
  if (p == "WHR") {
    colnames(res) <- c("eid", 
                       "event_dt.1", "event_dt.2", "event_dt.3",
                       "waist.1", "waist.2", "waist.3", 
                       "hip.1", "hip.2", "hip.3",
                       "data_provider.1", "data_provider.2", 
                       "data_provider.3",
                       "BMI.1", "BMI.2", "BMI.3")
    res <- pivot_longer(res, 
                        cols = -1, 
                        names_to = c(".value", "set"),
                        names_pattern = "(.+).(.+)", 
                        values_drop_na = T) %>%
      mutate(value = waist/hip)
  } else {
    # Match column names with GP data, main phenos have 2 repeats
    colnames(res) <- c("eid", 
                       "event_dt.1", "event_dt.2", "event_dt.3",
                       "value.1", "value.2", "value.3",
                       "data_provider.1", "data_provider.2", 
                       "data_provider.3",
                       "BMI.1", "BMI.2", "BMI.3")
    res <- pivot_longer(res, cols = -1, 
                        names_to = c(".value", "set"),
                        names_pattern = "(.+).(.+)", 
                        values_drop_na = T)
    
  } 
  
  res$data_provider <- paste0("UKBB", res$data_provider) 
  res$biomarker <- p
  res <- res[, c("eid", "data_provider", "event_dt",
                 "value", "BMI", "biomarker")]
  return (res)
}

# Function to get age at event from birth date
get_age <- function (long_df) {
  # Get date of birth from main phenotype file (month and year)
  tmp <- merge(long_df, pheno[, c("eid", "f.34.0.0", "f.52.0.0")],
               by = "eid")
  colnames(tmp) <- c(colnames(long_df), "year_of_birth", "month_of_birth")
  # Assign date of birth (first of the month) based on month and year of birth 
  tmp$dob <- as.Date(paste(tmp$year_of_birth, tmp$month_of_birth,
                           1, sep = "-"))
  
  # Calculate age at event
  tmp$event_dt <- as.Date(tmp$event_dt, "%Y-%m-%d")
  tmp$dob <- as.Date(tmp$dob, "%Y-%m-%d")
  # Remove inconsistencies (NA in event date, event date before date-of-birth, 
  # and event date after 2020)
  tmp <- tmp %>% filter(!is.na(event_dt) & event_dt > dob & 
                          event_dt <= as.Date("2020-10-01"))
  tmp$age_event <- interval(tmp$dob, tmp$event_dt) / years(1)
  
  res <- tmp[, c("eid", "data_provider", "event_dt", "age_event",
                 "value", "BMI", "biomarker")]
  return (res)
}

# Function to standardise column types 
standardise_coltypes <- function (long_df) {
  res <- long_df[, c("eid", "data_provider", 
                     "event_dt", "age_event", "value", "BMI",
                     "biomarker")]
  # Convert all columns to standard type
  res$eid <- factor(long_df$eid)
  res$data_provider <- factor(long_df$data_provider)
  res$event_dt <- as.Date(long_df$event_dt, "%Y-%m-%d")
  res$age_event <- as.numeric(long_df$age_event)
  res$BMI <- as.numeric(long_df$BMI)
  return (res)
}

# Function to get BMI-adjusted value
adj_for_bmi <- function (long_df) {
  complete_dat <- long_df[complete.cases(long_df), ]
  modeled_dat <- lm(value ~ BMI,
                    data = complete_dat)
  resids <- resid(modeled_dat)
  
  res <- complete_dat
  res$value_adjBMI <- resids
  
  # Add back in rows that were missing BMI
  add_rows <- which(is.na(long_df$BMI))
  add_rows <- long_df[add_rows, ]
  add_rows$value_adjBMI <- NA
  
  res <- bind_rows(res, add_rows)
  return (res)
}

all_adipo_main_dat <- lapply(PHENOTYPES, function (p) {
  print(paste0("Running phenotype: ", p))
  
  # Special case for WHR as it has to be calculated from WC and HC
  if (p == "WHR") {
    value_codes <- c("f.48.0.0", "f.48.1.0", "f.48.2.0",
                     "f.49.0.0", "f.49.1.0", "f.49.2.0")
    date_codes <- c("f.53.0.0", "f.53.1.0", "f.53.2.0")
    main_df <- pheno[, c("eid", date_codes, value_codes)]
  } else {
    # Get relevant columns for value, date, and assessment centre
    value_codes <- 
      colnames(pheno)[grep(PHENO_CODES[p],
                           colnames(pheno))]
    date_codes <- c("f.53.0.0", "f.53.1.0", "f.53.2.0")
    main_df <- pheno[, c("eid", date_codes, value_codes)]
  }
  # Convert main df to long format 
  main_df <- wide_to_long(main_df, p)
  # Get age at event from birth date
  main_df <- get_age(main_df)
  # Standardise columns
  res <- standardise_coltypes(main_df)
  
  # Get BMI-adjusted data where appropriate
  if (p %in% c("WC", "WHR")) {
    res <- adj_for_bmi(res)
  }
  
  return (res)
})
names(all_adipo_main_dat) <- PHENOTYPES

# Change coding on weight-change-1yr question
tmp <- all_adipo_main_dat[["Weight_change_1yr"]]
tmp$value <- ifelse(tmp$value == 0, "No change",
                    ifelse(tmp$value == 2, "Gain",
                           ifelse(tmp$value == 3, "Loss", "Unknown")))
tmp$value <- factor(tmp$value)
all_adipo_main_dat[["Weight_change_1yr"]] <- tmp

# Filter to longitudinal data except for the weight-change field ----

filtered <- lapply(PHENOTYPES, function (p) {
  res <- all_adipo_main_dat[[p]]
  res <- res %>% arrange(eid, event_dt)
  
  # Drop any NA values and only keep individuals with > 1 measurement
  cleaned <- res[complete.cases(res$age_event, res$value), ]
  keep_ids <- cleaned %>% group_by(eid) %>% count()

  if (p == "Weight_change_1yr") 
    keep_ids <- keep_ids$eid
  else 
    keep_ids <- keep_ids$eid[keep_ids$n > 1]
  
  cleaned <- subset(cleaned, cleaned$eid %in% keep_ids)
 
  return (cleaned)
})
names(filtered) <- PHENOTYPES

saveRDS(filtered, 
        "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/data/main_data_adipo_change.rds")

