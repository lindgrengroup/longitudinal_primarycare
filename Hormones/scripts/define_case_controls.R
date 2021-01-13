# Author: Samvida S. Venkatesh
# Date: 12/01/2021

library(tidyverse)
library(lubridate)

set.seed(120121)
CC_RATIO = 4 # case-control ratio 1:4

# Read data ----

# Read gp_clinical file
gp_clinical <- readRDS("/well/lindgren/UKBIOBANK/samvida/gp_clinical_annotated_sex_dob.rds")

# Keep demographic information
demo_info <- gp_clinical %>% distinct(eid, sex, dob, mean_UKBB_BMI)
demo_info$case <- F

# Relevant read codes for cases

# LH: plasma and serum
case_read2_codes <- c("443f.", "443e.", "4433")
case_read3_codes <- c("XaELa", "XM0lv", "X80Fu", "XE25I")

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
# a. DOB (+/- 6 months) & 
# b. UKBIOBANK BMI (+/- 1 unit of BMI)

############ WRITE A PIECE HERE TO FIGURE OUT HOW WIDE THE MATCHING CRITERIA #
# SHOULD BE BASED ON NUMBER OF MATCHES PER CASE ##############################

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
    control_demo <- subset(non_case_demo, 
                           non_case_demo$sex == case_sex &
                             non_case_demo$dob %within% case_dob_interval)
  } else {
    control_demo <- subset(non_case_demo, 
                           non_case_demo$sex == case_sex &
                             (non_case_demo$mean_UKBB_BMI >= (case_BMI*0.95) & 
                                non_case_demo$mean_UKBB_BMI <= (case_BMI*1.05)) &
                             (non_case_demo$dob %within% case_dob_interval))
  }
  matched_eids <- unique(control_demo$eid)
  n_matches <- length(matched_eids)
  
  # Sample from matched population for controls
  res <- sample(matched_eids, size = min(n_matches, CC_RATIO), replace = F)
  return (res)
}

control_eids_s3 <- lapply(case_eids, function (id) findMatchedControls(id) )
control_eids_s3 <- unique(unlist(control_eids_s3))

# Save list of IDs for downstream analyses
LH_eids <- list(case_eids,
                control_eids_s1,
                control_eids_s2,
                control_eids_s3)
saveRDS(LH_eids, "/well/lindgren/UKBIOBANK/samvida/hormone_ehr/LH/case_control_ids.rds")