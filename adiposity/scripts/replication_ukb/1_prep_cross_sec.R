# Author: Samvida S. Venkatesh
# Date: 18/08/22

library(tidyverse)

# Read data ----

PHENOTYPES <- c("BMI", "Weight")

# Main data
dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/main_data_without_gp_cross_sec.rds")[PHENOTYPES]

# GP data used for intercept-GWAS
gp_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/indiv_qcd_data.rds")[PHENOTYPES]

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, comment.char = "$",
                             stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

SEX_STRATA <- c("F", "M", "sex_comb")
# Add in data provider as covariate later if there are enough levels
COVARS <- c("age_event", "age_sq", "year_of_birth", "smoking_status")

# QC log file
qc_log <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/qc/cross_sectional_qc.txt"

# Wrangle data ----

# Ensure that there are no individuals with GP data in the no-GP file
dat <- lapply(PHENOTYPES, function (p) {
  res <- dat[[p]]
  remove_ids <- as.character(unique(gp_dat[[p]]$eid))
  
  res <- res %>% 
    mutate(eid = as.character(eid)) %>%
    filter(!eid %in% remove_ids)
  
  return (res)
})
names(dat) <- PHENOTYPES

# Get cross-sectional value
cross_sec_dat <- lapply(PHENOTYPES, function (p) {
  res <- dat[[p]] %>%
    group_by(eid) %>%
    # Get median value and distance of each observation to median
    mutate(medn_value = median(value),
           dist_to_medn = abs(value - median(value))) %>%
    # Arrange by age to get first observation that is closest to median
    group_by(eid) %>% arrange(age_event, .by_group = T) %>%
    slice(which.min(dist_to_medn))
  
  res <- res %>% 
    select(all_of(c("eid", "data_provider", "event_dt",
                    "age_event", "value", "biomarker")))
  
  return (res)
})
names(cross_sec_dat) <- PHENOTYPES

# Add covariates, adjust and RINT values ----

cross_sec_dat <- lapply(PHENOTYPES, function (p) {
  
  df <- cross_sec_dat[[p]]
  
  df <- left_join(df, 
                  general_covars[, c("eid", "sex", "year_of_birth", "smoking_status")],
                  by = "eid")
  
  # Only retain individuals that have covariate information 
  df <- df[complete.cases(df), ]
  sink(qc_log, append = T)
  cat(paste0("** PHENOTYPE ** ", p, "\n",
             "\t", "# males: ", sum(df$sex == "M"), "\n",
             "\t", "# females: ", sum(df$sex == "F"), "\n"))
  sink()
  
  # Calculate adjustments and transform (RINT) the value
  res <- lapply(SEX_STRATA, function (sx) {
    if (sx == "sex_comb") {
      sub_df <- df
      covars_adj <- c(COVARS, "sex")
    }
    else {
      sub_df <- df %>% filter(sex == sx)
      covars_adj <- COVARS
    }
    
    # Add covariate for squared age
    sub_df <- sub_df %>% mutate(age_sq = age_event^2)
    # Add data provider as a covariate if there is > 1 provider
    if (length(unique(sub_df$data_provider)) > 1) {
      covars_adj <- c(covars_adj, "data_provider")
    }
    
    # Formula for adjustment
    mod_formula <- 
      formula(paste0("value ~ ", paste(covars_adj, collapse = " + ")))
    # Get residuals
    mod_resid <- lm(mod_formula, data = sub_df)$residuals
    # RINT
    sub_df$rinted_value <- qnorm((rank(mod_resid) - 0.5) / sum(!is.na(mod_resid)))
    to_write <- data.frame(FID = sub_df$eid, 
                           IID = sub_df$eid,
                           adj_trait = sub_df$rinted_value)
    # Write results to table
    write.table(to_write,
                paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/GWAS/traits_for_gwas/", 
                       p, "_", sx, "_cross_sectional.txt"),
                sep = "\t", row.names = F, quote = F)
    
    return (sub_df)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(cross_sec_dat) <- PHENOTYPES

saveRDS(cross_sec_dat, 
        "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/data/split_sex_cross_sec.rds")
