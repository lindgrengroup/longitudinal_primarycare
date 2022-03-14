# Author: Samvida S. Venkatesh
# Date: 14/03/22

# Create new UKBB HES file with age at each event

# NEEDS AT LEAST 5 CORES (PREFERABLY 7 CORES)

# NPAR_CORES <- 7

library(lubridate)
library(eeptools)
# library(foreach)
# library(doParallel) 
library(dplyr)

# Read data ----

# Read HES file
hes_dat <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/HES/ukb_hesin_plus_birth_secdiag_secoper.tsv",
                      sep = "\t", header = T, comment.char = "~",
                      na.string = c("NA", "", "."),
                      stringsAsFactors = F)
ORIG_DAT_COLS <- colnames(hes_dat)

# Main UKBB phenotypes file
pheno <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "f.eid"

# Add info from main UKBB phenotypes file ----

## Add information on date-of-birth and sex ----

hes_dat <- merge(hes_dat, pheno[, c("f.eid", "f.31.0.0", 
                                    "f.34.0.0", "f.52.0.0")], 
                 by.x = "eid", by.y = "f.eid", all.x = T)
NCOL <- ncol(hes_dat)
colnames(hes_dat)[(NCOL-2):NCOL] <- c("sex", "year_of_birth", "month_of_birth")

# Assign date of birth (first of the month) based on month and year of birth 
hes_dat$dob <- as.Date(paste(hes_dat$year_of_birth, 
                             hes_dat$month_of_birth,
                             1, sep = "-"), "%Y-%m-%d")
# Remove records without DOB
hes_dat <- hes_dat %>% filter(!is.na(dob))

# Change sex coding to M, F
hes_dat$sex <- ifelse(hes_dat$sex == 0, "F", "M")

# Columns for which to calculate age ----

# admidate = date of admission to hospital
# disdate = date of discharge from hopsital 
# elecdate = date of decision to admit to hospital 
# epistart = episode start date
# epiend = episode end date 
# opdate = date of operation
# opdate_sec = date of operation (sec)

DATE_COLS_INTEREST <- c("admidate", "disdate",
                        "elecdate", 
                        "epistart", "epiend",
                        "opdate", "opdate_sec")
# Ensure that the columns are dates and remove inconsistencies

removeInconsistentDates <- function (date_x, dob_x) {
  # Convert to date
  res <- as.Date(date_x, "%Y-%m-%d")
  # Look for inconsistencies
  na_date <- which(is.na(res))
  date_before_dob <- which(res < dob_x)
  date_after_2020 <- which(res > as.Date("2020-10-01", "%Y-%m-%d"))
  # Replace inconsistencies with NA
  replace_with_na <- unique(c(na_date, date_before_dob, date_after_2020))
  res[replace_with_na] <- NA
  return (res)
}

hes_dat <- hes_dat %>% mutate(across(all_of(DATE_COLS_INTEREST), 
                                     ~ removeInconsistentDates(date_x = ., 
                                                               dob_x = dob)))

# Calculate age at each event ----

# Wrapper function for age-calc that ignores NAs
wrapAgeCalc <- function (date_x, dob_x) {
  to_calc <- which(!is.na(date_x))
  res <- rep(NA, length(date_x))
  # Calculate ages
  res[to_calc] <- age_calc(dob = dob_x[to_calc],
                           enddate = date_x[to_calc],
                           units = "years")
  return (res)
}

# Calculate age at each event
tmp <- hes_dat %>% mutate(across(all_of(DATE_COLS_INTEREST),
                                 ~ wrapAgeCalc(date_x = .,
                                               dob_x = dob),
                                 .names = "age_{.col}"))

# Get earliest date / youngest age in each row, across all types of dates
AGE_COLS_INTEREST <- paste0("age_", DATE_COLS_INTEREST)
tmp <- tmp %>% 
  rowwise() %>%
  mutate(age_record = min(c_across(all_of(AGE_COLS_INTEREST)), na.rm = T)) %>%
  # quality control age record by making sure it is > 0
  filter(age_record >= 0)

res <- tmp %>% select(all_of(c(ORIG_DAT_COLS, "dob", "age_record")))

write.table(res, 
            "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/hes_age_annotated.txt",
            sep = "\t", quote = F, row.names = F)
