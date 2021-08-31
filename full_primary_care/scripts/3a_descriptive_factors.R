# Author: Samvida S. Venkatesh
# Date: 27/07/21

library(tidyverse)
library(lubridate)

set.seed(020321)

# Read all biomarkers and covariate data ----

dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data_passed_popn_QC.rds")
PHENOTYPES <- names(dat)

general_covars <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/QCd_demographic_covariates.txt",
                             sep = "\t", header = T, comment.char = "$",
                             stringsAsFactors = F)

# Calculate trait-specific covariates ----

trait_covars <- lapply(PHENOTYPES, function (p) {
  res <- dat[[p]] %>% group_by(eid) %>% 
    arrange(eid, event_dt) %>%
    summarise(baseline_date = first(event_dt),
              baseline_age = first(age_event),
              FUyrs = interval(first(event_dt), last(event_dt)) / years(1),
              FU_n = n(),
              baseline_trait = first(value))
})
names(trait_covars) <- PHENOTYPES

# Construct demographic characteristics tables -----

demo_tables <- lapply(PHENOTYPES, function (p) {
  res <- merge(trait_covars[[p]], 
               general_covars[, c("eid", "sex", "ancestry", "height")],
               by = "eid")
  
  res <- res %>% group_by(sex) %>% 
    summarise(count = n(), 
              mean_FU_n = mean(FU_n), sd_FU_n = sd(FU_n),
              mean_FUyrs = mean(FUyrs), sd_FUyrs = sd(FUyrs),
              mean_bl_age = mean(baseline_age), 
              sd_bl_age = sd(baseline_age),
              median_bl_trait = median(baseline_trait), 
              iqr_bl_trait = paste(quantile(baseline_trait, 0.25), 
                                   quantile(baseline_trait, 0.75),
                                   sep = ", ")) %>%
    ungroup() %>% mutate(percent = count / sum(count))
  res$biomarker <- p
  return (res)
})
demo_tables <- bind_rows(demo_tables)

write.table(demo_tables, 
            paste0("results/primary_care_descriptive_factors.txt"),
            sep = "\t", quote = F, row.names = F)

