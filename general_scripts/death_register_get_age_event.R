# Author: Samvida S. Venkatesh
# Date: 26/07/22

# Create new UKBB death register file with age at each event

# NEED 7 CORES

library(lubridate)
library(eeptools)
library(dplyr)

# Read data ----

# Death dates file
death_dates <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/DATA/DEATH/death_21012022.txt",
                          sep = "\t", header = T, comment.char = "~",
                          na.string = c("NA", "", "."),
                          stringsAsFactors = F)
death_dates <- death_dates %>% 
  select(all_of(c("eid", "date_of_death"))) %>%
  distinct()

# Causes of death file
death_causes <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/DATA/DEATH/death_cause_21012022.txt",
                           sep = "\t", header = T, comment.char = "~",
                           na.string = c("NA", "", "."),
                           stringsAsFactors = F)
death_causes <- death_causes %>% 
  select(all_of(c("eid", "cause_icd10"))) %>%
  distinct()

# Combine the two by eid
death_dat <- full_join(death_dates, death_causes, by = "eid")

# Main UKBB phenotypes file
pheno <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "f.eid"

# Add info from main UKBB phenotypes file ----

## Add information on date-of-birth and sex ----

death_dat <- merge(death_dat, pheno[, c("f.eid", "f.31.0.0", 
                                        "f.34.0.0", "f.52.0.0")], 
                   by.x = "eid", by.y = "f.eid", all.x = T)
NCOL <- ncol(death_dat)
colnames(death_dat)[(NCOL-2):NCOL] <- c("sex", "year_of_birth", "month_of_birth")

# Assign date of birth (first of the month) based on month and year of birth 
death_dat$dob <- as.Date(paste(death_dat$year_of_birth, 
                               death_dat$month_of_birth,
                               1, sep = "-"), "%Y-%m-%d")
# Remove records without DOB
death_dat <- death_dat %>% filter(!is.na(dob))

# Change sex coding to M, F
death_dat$sex <- ifelse(death_dat$sex == 0, "F", "M")

# Columns for which to calculate age ----

# Ensure that the death date columns are dates and remove inconsistencies

removeInconsistentDates <- function (date_x, dob_x) {
  # Convert to date
  res <- as.Date(date_x, "%d/%m/%Y")
  # Look for inconsistencies
  na_date <- which(is.na(res))
  date_before_dob <- which(res < dob_x)
  date_after_2025 <- which(res > as.Date("2025-01-01", "%Y-%m-%d"))
  # Replace inconsistencies with NA
  replace_with_na <- unique(c(na_date, date_before_dob, date_after_2025))
  res[replace_with_na] <- NA
  return (res)
}

death_dat <- death_dat %>% 
  mutate(date_of_death = removeInconsistentDates(date_x = date_of_death, 
                                                 dob_x = dob))

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
tmp <- death_dat %>%
  mutate(age_at_death = wrapAgeCalc(date_x = date_of_death,
                                    dob_x = dob))

res <- tmp %>% 
  # quality control age record by making sure it is > 0 and not NA
  filter(!is.na(age_at_death) & age_at_death >= 0)

write.table(res, 
            "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/death_register_age_annotated.txt",
            sep = "\t", quote = F, row.names = F)
