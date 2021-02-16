# Author: Samvida S. Venkatesh
# Date: 13/02/21

library(tidyverse)
library(lubridate)

# Read data ----

# Read full primary care data
gp_clinical <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/gp_clinical_annotated.txt",
                          sep = "\t", header = T, comment.char = "$")
gp_clinical$event_dt <- as.Date(gp_clinical$event_dt, "%Y-%m-%d")
gp_clinical$dob <- as.Date(gp_clinical$dob, "%Y-%m-%d")

# Read full UKB phenotypes data
pheno <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "f.eid"

# Read pregnancy codes in primary care
pregnancy_codes <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/pregnancy_codes.txt",
                              sep = "\t", header = T, comment.char = "$")

# It is difficult to get pregnancy periods with much precision, so flag data from
# +/- 1 year around each pregnancy code

# Primary care pregnancy events ----

preg_EIDDATES_gp <- gp_clinical[gp_clinical$read_2 %in% 
                                  pregnancy_codes$read_code |
                                  gp_clinical$read_3 %in% 
                                  pregnancy_codes$read_code, c("eid", "event_dt")]

preg_EIDDATES_gp <- preg_EIDDATES_gp %>%
  mutate(start = event_dt - years(1), end = event_dt + years(1)) 

# UKB pregnancy records ----

EVENT_DATE_CODES <- c("f.53.0.0", "f.53.1.0", "f.53.2.0")
CURRENT_PREG_CODES <- c("f.3140.0.0", "f.3140.1.0", "f.3140.2.0")
PREG_AGE_CODES <- c("f.3872.0.0", "f.3872.1.0", "f.3872.2.0", # primiparous
                    "f.2754.0.0", "f.2754.1.0", "f.2754.2.0", # first birth
                    "f.2764.0.0", "f.2764.1.0", "f.2764.2.0") # last birth

## First get "current" pregnancy dates recorded in UKB ----
preg_EIDDATES_ukbb <- pheno[, c("f.eid", 
                                EVENT_DATE_CODES, 
                                CURRENT_PREG_CODES)]
colnames(preg_EIDDATES_ukbb) <- c("eid", 
                                  "event_dt.1", "event_dt.2", "event_dt.3",
                                  "current_preg.1", "current_preg.2", 
                                  "current_preg.3")
preg_EIDDATES_ukbb <- pivot_longer(preg_EIDDATES_ukbb, cols = -1, 
                                   names_to = c(".value", "set"),
                                   names_pattern = "(.+).(.+)", 
                                   values_drop_na = T)
# Ensure event date column is a date
preg_EIDDATES_ukbb$event_dt <- as.Date(preg_EIDDATES_ukbb$event_dt, "%Y-%m-%d")
# Subset pregnancies
preg_EIDDATES_ukbb <- subset(preg_EIDDATES_ukbb, 
                             preg_EIDDATES_ukbb$current_preg == 1)
preg_EIDDATES_ukbb <- preg_EIDDATES_ukbb[, c("eid", "event_dt")] %>%
  mutate(start = event_dt - years(1), end = event_dt + years(1)) 

## Then use recorded age at childbirth data ----

preg_EIDAGES_ukbb <- pheno[, c("f.eid", PREG_AGE_CODES)]
colnames(preg_EIDAGES_ukbb) <- c("eid", 
                                  "age_preg1.1", "age_preg1.2", "age_preg1.2",  
                                  "age_pregf.1", "age_pregf.2", "age_pregf.3", 
                                  "age_pregl.1", "age_pregl.2", "age_pregl.3")
preg_EIDAGES_ukbb <- pivot_longer(preg_EIDAGES_ukbb, cols = -1, 
                                   names_to = c(".value", "set"),
                                   names_pattern = "(.+).(.+)", 
                                  values_drop_na = T)
# Calculate median age from the three visits (ex. if recalled incorrectly) 
# NOT ACROSS DIFFERENT CHILDREN
preg_EIDAGES_ukbb <- preg_EIDAGES_ukbb %>% group_by(eid) %>%
  summarise(age_preg1 = median(age_preg1, na.rm = T),
            age_pregf = median(age_pregf, na.rm = T),
            age_pregl = median(age_pregl, na.rm = T))
# Pivot into separate pregnancies per individual
preg_EIDAGES_ukbb <- pivot_longer(preg_EIDAGES_ukbb, cols = -1, 
                               names_to = "preg_type",
                               values_to = "preg_age", 
                               values_drop_na = T)
# Retain the right ages (not negative numbers which indicate missingness)
preg_EIDAGES_ukbb <- subset(preg_EIDAGES_ukbb, preg_EIDAGES_ukbb$preg_age > 0)

# Calculate pregnancy dates from ages at childbirth
EIDS_DOB <- gp_clinical %>% distinct(eid, dob)

preg_EIDAGES_ukbb$dob <- EIDS_DOB$dob[match(preg_EIDAGES_ukbb$eid, 
                                            EIDS_DOB$eid)]
preg_EIDAGES_ukbb$dob <- as.Date(preg_EIDAGES_ukbb$dob, "%Y-%m-%d")

preg_EIDAGES_ukbb <- preg_EIDAGES_ukbb %>% 
  mutate(start = dob + years(as.integer(preg_age - 1)), 
         end = dob + years(as.integer(preg_age + 1)))

# Combine UKB and GP data ----

preg_EIDDATES <- bind_rows(preg_EIDDATES_gp[, c("eid", "start", "end")],
                           preg_EIDDATES_ukbb[, c("eid", "start", "end")],
                           preg_EIDAGES_ukbb[, c("eid", "start", "end")])

# Flatten intervals by combining adjacent sequential intervals
# CODE ADAPTED FROM: 
# https://stackoverflow.com/questions/28938147/how-to-flatten-merge-overlapping-time-periods

# Ensure ordered by pregnancy record date
preg_EIDDATES <- preg_EIDDATES %>% 
  arrange(eid, start) %>%
  group_by(eid) %>%
  # Calculate index for each pregnancy
  mutate(ind = c(0, cumsum(as.numeric(lead(start)) >
                             cummax(as.numeric(end)))[-n()])) %>%
  # Interval within each pregnancy and individual
  group_by(eid, ind) %>%
  summarise(start = first(start), end = last(end)) 

write.table(preg_EIDDATES, "/well/lindgren/UKBIOBANK/samvida/general_resources/pregnancy_records.txt",
            sep = "\t", row.names = F, quote = F)
