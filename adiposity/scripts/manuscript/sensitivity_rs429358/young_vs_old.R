# Author: Samvida S. Venkatesh
# Date: 05/10/21

library(lme4)
library(tidyverse)
theme_set(theme_bw())

# Read files ----

mainpath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data"
gpdat_path <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity"
gen_resources_path <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources"
outpath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sensitivity_rs429358"

# Parse in phenotype argument
# args = commandArgs(trailingOnly=TRUE)
PHENO <- "Weight"
SEX_STRATA <- c("F", "M", "sex_comb")

dat <- readRDS(paste0(mainpath, "/indiv_qcd_data.rds"))[[PHENO]]
dat$eid <- as.character(dat$eid)
covars <- readRDS(paste0(mainpath, "/covariates.rds"))[[PHENO]]
covars$eid <- as.character(covars$eid)

general_covars <- read.table(paste0(gen_resources_path, "/220504_QCd_demographic_covariates.txt"),
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

# Genotypes / dosages at rs429358
apoe_dosages <- read.table(paste0(gpdat_path, "/sample_variant_counts/rs429358_dosages.txt"),
                           sep = " ", header = T, stringsAsFactors = F)
# Remove first row, which contains info on type of column and columns 
# 2, 3, 4 (ID repeat, missingness, sex)
apoe_dosages <- apoe_dosages[-1, c(1, 5)]
colnames(apoe_dosages) <- c("eid", "dosage")
apoe_dosages$eid <- as.character(apoe_dosages$eid)

# Wrangle data ----

# Get time from baseline measurement and add in covariates
model_dat <- left_join(dat, 
                       covars[, c("eid", "baseline_age", "age_sq")], 
                       by = "eid")
model_dat <- left_join(model_dat, 
                       general_covars[, c("eid", "sex", "year_of_birth", paste0("PC", 1:21))],
                       by = "eid")
model_dat <- model_dat %>% 
  mutate(t = age_event - baseline_age)

# Add in APOE genotype dosage
model_dat <- model_dat %>%
  left_join(apoe_dosages, by = "eid") %>%
  mutate(dosage = as.numeric(dosage))

model_dat <- model_dat[complete.cases(model_dat), ]

# Run models ----

getSlope <- function (dat) {
  lmod <- lmer(value ~ t + (t | eid), 
               data = dat, REML = F)
  res <- coef(lmod)$eid
  # We only want the columns with both fixed and random effects 
  # (intercept and "t")
  # But we want to get the full effect (fixed + random) for it 
  # Return dataframe to write
  res$eid <- rownames(res)
  res <- res[, c("eid", "t")]
  colnames(res) <- c("eid", "weight_slope")
  return (res)
}

plotDist <- function (dat) {
  means_dat <- dat %>% group_by(type) %>%
    summarise(weight_slope = mean(weight_slope))
  
  resplot <- ggplot(dat, aes(x = weight_slope,
                             fill = type, colour = type)) +
    geom_point(data = means_dat, 
               aes(x = weight_slope, y = 0, colour = type)) +
    geom_density(alpha = 0.3) +
    scale_fill_brewer(palette = "Dark2") +
    scale_colour_brewer(palette = "Dark2")
  return (resplot)
}

full_models <- lapply(SEX_STRATA, function (sx) {
  by_geno <- lapply(c("all", "only_TT"), function (gty) {
    sub_dat <- model_dat
    if (sx != "sex_comb") {
      sub_dat <- sub_dat %>% filter(sex == sx)
    }
    if (gty == "only_TT") {
      sub_dat <- sub_dat %>% filter(dosage == 0)
    }
    
    # Split to "young" (values between ages 30-50) and "old" (valus between ages 60-80)
    young_dat <- sub_dat %>% filter(age_event >= 30 & age_event <= 50)
    old_dat <- sub_dat %>% filter(age_event >= 60 & age_event <= 80)
    
    young_slopes <- getSlope(young_dat)
    young_slopes$type <- "ages [30, 50]"
    old_slopes <- getSlope(old_dat)
    old_slopes$type <- "ages [60, 80]"
    
    dist_plot <- plotDist(bind_rows(young_slopes, old_slopes)) +
      labs(title = paste0(gty, "_", sx))
    ggsave(paste0(outpath, "/plot_young_vs_old_weight_slopes_", gty, "_", sx, ".png"),
           dist_plot, width = 7.5, height = 7.5, units = "in")
    
  })
})
