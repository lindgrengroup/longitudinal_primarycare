# Author: Samvida S. Venkatesh
# Date: 12/01/2021

library(tidyverse)
library(lubridate)

set.seed(120121)
CC_RATIO = 4 # case-control ratio 1:4

# Read data ----

# Read gp_clinical file
gp_clinical <- readRDS("/well/lindgren/UKBIOBANK/samvida/gp_clinical_annotated.rds")

# Keep demographic information
demo_info <- gp_clinical %>% distinct(eid, sex, dob, mean_UKBB_BMI)
demo_info$case <- F

# Relevant read codes for cases

# GROUP 1: LUTEINISING HORMONE AND FOLLICLE STIMULATING HORMONE
case_read2_codes <- c("46J6.", "46J7.", "4434", "4433", "443f.", "443i.",
                      "44342", "44332", "4Q23.", "443e.", "443h.")
case_read3_codes <- c("4433", "4434", "46J6.", "46J7.", "X7722", "X80Fa",
                      "X80Fu", "XaELa", "XaELZ", "XaX4U", "XaX4V", "XE25I",
                      "XE25J", "XM0lv", "XM0lx")

SEXES <- c("F", "M")

# Calculate case IDs ----

case_eids <- subset(gp_clinical, gp_clinical$read_2 %in% case_read2_codes | 
                      gp_clinical$read_3 %in% case_read3_codes)$eid
case_eids <- unique(case_eids)
demo_info$case[demo_info$eid %in% case_eids] <- T

demo_info <- demo_info %>% group_by(sex) %>% 
  mutate(NCASES = sum(case))

# Calculate control IDs ----

# Strategy 1: Random cohort

control_eids_s1 <- demo_info %>% group_by(sex) %>% 
  sample_n(size = NCASES * CC_RATIO, replace = F)
control_eids_s1 <- control_eids_s1$eid

# Strategy 2: Random non-case cohort

non_case_demo <- subset(demo_info, !demo_info$case)
control_eids_s2 <- non_case_demo %>% group_by(sex) %>% 
  sample_n(size = NCASES * CC_RATIO, replace = F)
control_eids_s2 <- control_eids_s2$eid

# Strategy 3: Matched controls
# a. DOB (+/- 1 year) & 
# b. UKBIOBANK BMI (+/- 5% of BMI)

findMatchedControls <- function (case_eid) {
  # Get case demographics
  case_demo <- subset(demo_info, demo_info$eid == case_eid)
  case_sex <- case_demo$sex
  case_BMI <- case_demo$mean_UKBB_BMI
  case_dob <- case_demo$dob
  case_dob_interval <- interval(case_dob %m-% months(12), case_dob %m+% months(12))
  
  # Subset primary care data with matching demographics
  if (is.na(case_BMI)) {
    # if we don't know case BMI
    control_demo <- subset(matched_controls_demo, 
                           matched_controls_demo$sex == case_sex &
                             matched_controls_demo$dob %within% case_dob_interval)
  } else {
    control_demo <- subset(matched_controls_demo, 
                           matched_controls_demo$sex == case_sex &
                             (matched_controls_demo$mean_UKBB_BMI >= (case_BMI*0.95) & 
                                matched_controls_demo$mean_UKBB_BMI <= (case_BMI*1.05)) &
                             (matched_controls_demo$dob %within% case_dob_interval))
  }
  matched_eids <- unique(control_demo$eid)
  
  # if not enough matches could be found based on both criteria then 
  # randomly sample from the right age group
  if (length(matched_eids) == 0) {
    control_demo <- subset(matched_controls_demo, 
                           matched_controls_demo$sex == case_sex &
                             matched_controls_demo$dob %within% case_dob_interval)
  }
  
  matched_eids <- unique(control_demo$eid)
  
  # Sample from matched population for controls
  res <- sample(matched_eids, 1)
  
  # Update matched control database to remove the chosen control
  matched_controls_demo <- subset(matched_controls_demo, 
                                  matched_controls_demo$eid != res)
  return (res)
}

matched_controls_demo <- non_case_demo
control_eids_s3 <- lapply(case_eids, function (id) findMatchedControls(id) )
control_eids_s3 <- unlist(control_eids_s3)

# Save list of IDs for downstream analyses
eids <- list(case_eids,
             control_eids_s1,
             control_eids_s2,
             control_eids_s3)
names(eids) <- c("case_eids", "control_eids_s1", 
                 "control_eids_s2", "control_eids_s3")

saveRDS(eids, "/well/lindgren/UKBIOBANK/samvida/hormone_ehr/LH_FSH/case_control_ids.rds")
