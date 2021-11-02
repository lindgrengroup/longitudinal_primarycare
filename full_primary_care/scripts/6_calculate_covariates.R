# Author: Samvida S. Venkatesh
# Date: 04/10/21

library(tidyverse)
library(lubridate)

# Read files ----

qcd_dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/indiv_qcd_data.rds")
general_covars <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/QCd_demographic_covariates.txt",
                             sep = "\t", header = T, comment.char = "$",
                             stringsAsFactors = F)
PHENOTYPES <- names(qcd_dat)

# QC log file
qc_log <- "/well/lindgren/UKBIOBANK/samvida/full_primary_care/qc/covariate_QC.txt"

# Design covariate files for each adiposity trait ----

NPCs = 21

QCd_covars <- lapply(PHENOTYPES, function (p) {
  df <- qcd_dat[[p]]
  # Calculate baseline and follow-up years covariates
  calc_covars <- df %>% group_by(eid) %>% 
    arrange(eid, event_dt) %>%
    summarise(baseline_date = first(event_dt),
              baseline_age = first(age_event),
              age_sq = baseline_age^2,
              FUyrs = interval(first(event_dt), last(event_dt)) / years(1),
              FU_n = n(),
              baseline_trait = first(value))
  
  # Merge with previously calculated covariates (sex, ancestry, etc.)
  cleaned <- merge(calc_covars, general_covars, by = "eid")
  # Set missing and inconsistent ancestry to "other"
  # Remove individuals missing any covariate
  res <- cleaned[, c("eid", "sex", "ancestry",
                     "baseline_age", "age_sq", 
                     "height", "baseline_trait",
                     "FUyrs", "FU_n", paste0("PC", 1:NPCs))]
  res <- res[complete.cases(res), ]
  sink(qc_log, append = T)
  cat(paste0("** PHENOTYPE **", p, "\n",  
             "**FILTER** EXCLUDED, Missing covariates: ", 
             dim(cleaned)[1] - dim(res)[1], "\n"))
  sink()
  
  res$eid <- as.character(res$eid)
  return(res)
})
names(QCd_covars) <- PHENOTYPES

# Save
saveRDS(QCd_covars, "/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/covariates.rds")
