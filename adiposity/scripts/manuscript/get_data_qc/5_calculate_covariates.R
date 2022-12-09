# Author: Samvida S. Venkatesh
# Date: 04/10/21

library(tidyverse)
library(lubridate)

mainpath <- "" # REDACTED

# Read files ----

qcd_dat <- readRDS(paste0(mainpath, "/data/indiv_qcd_data.rds"))
PHENOTYPES <- names(qcd_dat)

# QC log file
qc_log <- paste0(mainpath, "/qc/covariate_QC.txt")

# Design covariate files for each trait ----

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
  
  # Set missing and inconsistent ancestry to "other"
  # Remove individuals missing any covariate
  res <- calc_covars[, c("eid", "baseline_age", "age_sq", 
                     "baseline_trait",
                     "FUyrs", "FU_n")]
  res <- res[complete.cases(res), ]
  sink(qc_log, append = T)
  cat(paste0("** PHENOTYPE **", p, "\n",  
             "**FILTER** EXCLUDED, Missing covariates: ", 
             dim(calc_covars)[1] - dim(res)[1], "\n"))
  sink()
  
  res$eid <- as.character(res$eid)
  return(res)
})
names(QCd_covars) <- PHENOTYPES

# Save
saveRDS(QCd_covars, psate0(mainpath, "/data/covariates.rds"))
