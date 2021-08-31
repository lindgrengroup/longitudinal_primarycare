# Author: Samvida S. Venkatesh
# Date: 31/08/21

library(tidyverse)

# Read biomarker data ----

dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data_passed_popn_QC.rds")
PHENOTYPES <- names(dat)

general_covars <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/QCd_demographic_covariates.txt",
                             sep = "\t", header = T, comment.char = "$",
                             stringsAsFactors = F)

# Calculate sex-heterogeneity in baseline trait value ----

sex_het <- lapply(PHENOTYPES, function (p) {
  # Get baseline trait value
  df <- dat[[p]] %>% group_by(eid) %>% 
    arrange(eid, event_dt) %>%
    summarise(baseline_trait = first(value))
  # Add sex info
  df <- merge(df, general_covars[, c("eid", "sex")],
               by = "eid")
  
  # Calculate heterogeneity with two-sample t-test
  het_res <- t.test(baseline_trait ~ sex, data = df, 
                    alternative = "two.sided", mu = 0,
                    paired = F, var.equal = F, conf.level = 0.05)
  # Create results df
  res <- df %>% group_by(sex) %>% 
    summarise(mean_bl_trait = mean(baseline_trait),
              sd_bl_trait = sd(baseline_trait),
              n = n()) %>%
    mutate(biomarker = p)
  res$pval_het <- rep(het_res$p.value, 2)
  
  return (res)
  
})
sex_het <- bind_rows(sex_het)

write.table(sex_het, 
            "/well/lindgren/UKBIOBANK/samvida/full_primary_care/results/primary_care_baseline_sexhet.txt",
            sep = "\t", quote = F, row.names = F)

