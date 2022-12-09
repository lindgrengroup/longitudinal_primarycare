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
pheno$eid <- as.character(pheno$eid)

# Subset to non-WB ancestry IDs
wb_ids <- as.character(read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eids_white_british.txt",
                     header = F, sep = "\t", stringsAsFactors = F)$V1)

WTCHG_CODES <- colnames(pheno)[grep("f.2306.*", colnames(pheno))]
DATE_CODES <- c("f.53.0.0", "f.53.1.0", "f.53.2.0")
CTR_CODES <- c("f.54.0.0", "f.54.1.0", "f.54.2.0")
BMI_CODES <- c("f.21001.0.0", "f.21001.1.0", "f.21001.2.0")

# Subset to non-WB ancestry and get self-rep weight-change data ----
# Also retain BMI in order to adjust for BMI later

non_wb_subset <- pheno %>%
  filter(!eid %in% wb_ids)

# Get relevant columns for value, date, assessment centre, and BMI
cleaned_res <- non_wb_subset[, c("eid", DATE_CODES, WTCHG_CODES,
                                 CTR_CODES, BMI_CODES)]
# Convert to long format 
# Match column names with GP data, main phenos have 2 repeats
colnames(cleaned_res) <- c("eid", 
                   "event_dt.1", "event_dt.2", "event_dt.3",
                   "selfrep_wtchg.1", "selfrep_wtchg.2", "selfrep_wtchg.3",
                   "data_provider.1", "data_provider.2", 
                   "data_provider.3",
                   "BMI.1", "BMI.2", "BMI.3")
cleaned_res <- pivot_longer(cleaned_res, cols = -1, 
                    names_to = c(".value", "set"),
                    names_pattern = "(.+).(.+)", 
                    values_drop_na = T)
cleaned_res$data_provider <- paste0("UKBB", cleaned_res$data_provider) 
cleaned_res <- cleaned_res[, c("eid", "data_provider", "event_dt",
               "selfrep_wtchg", "BMI")]

# Get age at event from birth date ----

# Get date of birth from main phenotype file (month and year)
tmp <- merge(cleaned_res, non_wb_subset[, c("eid", "f.34.0.0", "f.52.0.0")],
             by = "eid")
colnames(tmp) <- c(colnames(cleaned_res), "year_of_birth", "month_of_birth")
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

cleaned_res <- tmp[, c("eid", "data_provider", "event_dt", "age_event",
               "selfrep_wtchg", "BMI")]

# Standardise columns ----

cleaned_res <- cleaned_res[, c("eid", "data_provider", 
                   "event_dt", "age_event", "selfrep_wtchg", "BMI")]
# Convert all columns to standard type
cleaned_res$eid <- factor(cleaned_res$eid)
cleaned_res$data_provider <- factor(cleaned_res$data_provider)
cleaned_res$event_dt <- as.Date(cleaned_res$event_dt, "%Y-%m-%d")
cleaned_res$age_event <- as.numeric(cleaned_res$age_event)
cleaned_res$BMI <- as.numeric(cleaned_res$BMI)

# Change coding on weight-change
cleaned_res <- cleaned_res %>% 
  mutate(selfrep_wtchg = ifelse(selfrep_wtchg == 0, "No change",
                                ifelse(selfrep_wtchg == 2, "Gain",
                                       ifelse(selfrep_wtchg == 3, "Loss", "Unknown"))))
write.table(cleaned_res, 
        "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/data/selfrep_wtchg_non_wb.txt",
        sep = "\t", row.names = F, quote = F)

