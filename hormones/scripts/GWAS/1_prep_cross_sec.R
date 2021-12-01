# Author: Samvida S. Venkatesh
# Date: 16/11/2021

library(tidyverse)

# Read data ----

HORMONES <- read.table("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/hormone_list.txt")$V1
dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/data_popn_qcd_no_longit_filter.rds")[HORMONES]
general_covars <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/QCd_demographic_covariates.txt",
                             sep = "\t", header = T, comment.char = "$",
                             stringsAsFactors = F)

SEX_STRATA <- c("F", "M", "sex_comb")
# Add in data provider as covariate later if there are enough levels
COVARS <- c("age_event", "age_sq")

# QC log file
qc_log <- "/well/lindgren/UKBIOBANK/samvida/hormone_ehr/qc/cross_sectional_qc.txt"

# Get cross-sectional value for each individual ----

# Define cross-sectional value as the first observed value closest to
# an individual's median and get age at event for this observed value

cross_sec_dat <- lapply(HORMONES, function (hr) {
  res <- dat[[hr]] %>% group_by(eid) %>%
    # Get median value and distance of each observation to median
    mutate(medn_value = median(value),
           dist_to_medn = abs(value - median(value))) %>%
    # Arrange by age to get first observation that is closest to median
    group_by(eid) %>% arrange(age_event, .by_group = T) %>%
    slice(which.min(dist_to_medn)) 
  return (res)
})
names(cross_sec_dat) <- HORMONES

# Add sex as covariate, adjust and RINT values ----

cross_sec_dat <- lapply(HORMONES, function (hr) {
  
  df <- cross_sec_dat[[hr]]
  df$sex <- general_covars$sex[match(df$eid, general_covars$eid)]
  # Only retain individuals that have sex information 
  df <- df[complete.cases(df), ]
  sink(qc_log, append = T)
  cat(paste0("** HORMONE ** ", hr, "\n",
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
    # Write formatted data to file for GWAS if there are enough, i.e. > 1000
    if (nrow(sub_df) > 1000) {
      to_write <- data.frame(FID = sub_df$eid, 
                             IID = sub_df$eid,
                             adj_trait = sub_df$rinted_value)
      # Write results to table
      write.table(to_write,
                  paste0("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS/traits_for_gwas/", 
                         hr, "_", sx, "_cross_sectional.txt"),
                  sep = "\t", row.names = F, quote = F)
    }
    return (sub_df)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(cross_sec_dat) <- HORMONES

saveRDS(cross_sec_dat, 
        "/well/lindgren/UKBIOBANK/samvida/hormone_ehr/data/split_sex_cross_sec.rds")