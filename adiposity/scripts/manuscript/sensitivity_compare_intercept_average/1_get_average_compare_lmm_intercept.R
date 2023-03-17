# Author: Samvida S. Venkatesh
# Date: 05/10/21

library(lme4)
library(tidyverse)
theme_set(theme_bw())

# Read files ----

mainpath <- "" # REDACTED
gen_resources_path <- "" # REDACTED
lmm_mods_path <- "" # REDACTED
results_path <- "" # REDACTED

PHENO <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

# Add sex as covariate for sex-combined analyses
MOD_COVARS <- c("baseline_age", "age_sq")
# Add data provider as covariate if there is more than one data provider
ADD_COVARS <- c("year_of_birth")

dat <- lapply(PHENO, function (ph) {
  res <- readRDS(paste0(mainpath, "/indiv_qcd_data.rds"))[[ph]]
  res$eid <- as.character(res$eid)
  return (res)
})
names(dat) <- PHENO

covars <- lapply(PHENO, function (ph) {
  res <- readRDS(paste0(mainpath, "/covariates.rds"))[[ph]]
  res$eid <- as.character(res$eid)
  return (res)
})
names(covars) <- PHENO

general_covars <- read.table(paste0(gen_resources_path, "/220504_QCd_demographic_covariates.txt"),
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

blups <- lapply(PHENO, function (ph) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    res <- read.table(paste0(lmm_mods_path, "/", ph, "_", sx, "_blups_full_model.txt"),
                      sep = "\t", header = T, stringsAsFactors = F)
    res$eid <- as.character(res$eid)
    colnames(res) <- c("eid", "b0", "b1")
    return (res)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(blups) <- PHENO

# Wrangle data ----

# Get average per individual, time from baseline measurement and add in covariates
avg_dat <- lapply(PHENO, function (ph) {
  res <- dat[[ph]] %>%
    summarise(mean_trait = mean(value))

  write.table(res, 
              paste0(results_path, "/", ph, "_average_value.txt"),
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  return (res)
})
names(avg_dat) <- PHENO

# Plot comparison between average BMI and BLUP from model ----

lapply(PHENO, function (ph) {
  lapply(SEX_STRATA, function (sx) {
    for_comp <- blups[[ph]][[sx]] %>%
      left_join(avg_dat[[ph]], by = "eid")
    
    corr_forplot <- round(cor(for_comp$mean_trait, for_comp$b0), 3)
    
    comp_plot <- ggplot(for_comp, aes(x = mean_trait, y = b0)) +
      geom_point(size = 0.3) +
      geom_smooth(method = lm) +
      labs(title = paste0("R2: ", corr_forplot))
    
    ggsave(paste0(results_path, "/plot_avg_", ph, "_", sx, "_vs_lmm_b0.png"),
           comp_plot, height = 7, width = 7, units = "in")
  })
})

