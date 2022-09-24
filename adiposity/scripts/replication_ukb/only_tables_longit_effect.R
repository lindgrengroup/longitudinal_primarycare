# Author: Samvida S. Venkatesh
# Date: 31/08/22

library(tidyverse)
library(broom)
library(foreign)
library(MASS)
theme_set(theme_bw())

# Read in arguments ----

resdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/replication_results/"

# Read data ----

main_longit_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/data/main_data_adipo_change.rds")
QUANT_PHENOS <- c("BMI", "Weight", "WC", "WHR",
                  "WCadjBMI", "WHRadjBMI")
CAT_PHENOS <- "Weight_change_1yr"

gp_ids <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/data/all_ids_in_discovery_gp_ukb_dat.txt",
                     sep = "\t", header = F, stringsAsFactors = F)$V1
gp_ids <- as.character(gp_ids)

vars_to_replicate <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/data/adipo_change_snps_replicate.txt",
                                sep = "\t", header = T, stringsAsFactors = F)
VARIDS <- vars_to_replicate$SNP

# Genotypes / dosages at variants of interest
var_dosages <- lapply(VARIDS, function (varid) {
  res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sample_variant_counts/",
                           varid, "_dosages.txt"),
                    sep = " ", header = T, stringsAsFactors = F)
  # Remove first row, which contains info on type of column and columns 
  # 2, 3, 4 (ID repeat, missingness, sex)
  res <- res[-1, c(1, 5)]
  colnames(res) <- c("eid", varid)
  return (res)
})
names(var_dosages) <- VARIDS

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

MOD_COVARS <- c("sex", "baseline_age", "age_sq", "FUyrs",
                "year_of_birth", "data_provider", 
                paste0("PC", 1:21))
general_covars <- general_covars %>% 
  dplyr::select(any_of(c("eid", MOD_COVARS))) %>%
  mutate(sex = factor(sex), year_of_birth = as.numeric(year_of_birth))

SEX_STRATA <- c("F", "M", "sex_comb")

# Wrangle data ----

# Convert dosages to 0/1/2 genotypes based on threshold
DOSAGE_THRESHOLD_0 <- 0.5
DOSAGE_THRESHOLD_2 <- 1.5
var_dosages_hardcall <- lapply(VARIDS, function (varid) {
  res <- var_dosages[[varid]] %>% 
    mutate(genotype = ifelse(!!as.symbol(varid) < DOSAGE_THRESHOLD_0, "0",
                             ifelse(!!as.symbol(varid) > DOSAGE_THRESHOLD_2, 
                                    "2", "1")),
           eid = as.character(eid)) %>%
    dplyr::select(all_of(c("eid", "genotype")))
  return (res)
})
names(var_dosages_hardcall) <- VARIDS

# Remove IDs with GP data from main longitudinal data for replication
# Tally visit number and changes from first visit needed for model
# Add in model covariates

dat_for_replication <- lapply(c(QUANT_PHENOS, CAT_PHENOS), function (p) {
  dat <- main_longit_dat[[p]] %>% 
    filter(!eid %in% gp_ids) %>%
    group_by(eid) %>% 
    mutate(visit = seq(n()),
           baseline_age = first(age_event),
           age_sq = baseline_age^2,
           age_event_sq = age_event^2,
           FUyrs = age_event - first(age_event)) %>%
    left_join(general_covars)
  
  if (p %in% QUANT_PHENOS) 
    dat <- dat %>% mutate(value_diff = value - first(value))
  else if (p %in% CAT_PHENOS)
    dat <- dat %>% mutate(value = factor(value, 
                                         levels = c("Loss", "No change", "Gain"),
                                         ordered = T))
  
  dat <- dat[complete.cases(dat), ]
  return (dat)
})
names(dat_for_replication) <- c(QUANT_PHENOS, CAT_PHENOS)

# Function to add variant dosage
addVarDosage <- function (dat, varid) {
  res <- dat
  var_dat <- var_dosages[[varid]]
  res$dosage <- var_dat[match(res$eid, var_dat$eid), varid]
  res$dosage <- as.numeric(res$dosage)
  res <- res[!is.na(res$dosage), ]
  return (res)
}

# Function to add hard-called genotype
addGenoGroup <- function (dat, varid) {
  res <- dat
  var_dat <- var_dosages_hardcall[[varid]]
  res$genotype <- var_dat$genotype[match(res$eid, var_dat$eid)]
  res$genotype <- factor(res$genotype, levels = c("0", "1", "2"))
  res <- res[!is.na(res$genotype), ]
  return (res)
}

# Testing functions ----

quantTest <- function (dat, ss) {
  covars_include <- MOD_COVARS
  if (ss != "sex_comb") {
    dat <- dat %>% filter(sex == ss)
    covars_include <- covars_include[-which(covars_include == "sex")]
  }
  if (length(unique(dat$data_provider)) < 2)
    covars_include <- covars_include[-which(covars_include == "data_provider")]
  
  mod_formula <- paste0("value_diff ~ dosage + ",
                        paste0(covars_include, collapse = " + "))
  
  modeled_dat <- lm(formula(mod_formula), data = dat)
  print_res <- tidy(modeled_dat) %>%
    filter(term == "dosage") %>%
    rename(beta = estimate, se = std.error, tstat = statistic, pval = p.value) %>%
    dplyr::select(all_of(c("beta", "se", "tstat", "pval"))) %>%
    mutate(sample_size = nrow(dat))
  
  return (print_res)
}

catTest <- function (dat, ss) {
  # Get the correct covariates for adjustment
  covars_include <- c("sex", "age_event", "age_event_sq",
                      "year_of_birth", "data_provider", 
                      paste0("PC", 1:21))
  if (ss != "sex_comb") {
    dat <- dat %>% filter(sex == ss)
    covars_include <- covars_include[-which(covars_include == "sex")]
  }
  if (length(unique(dat$data_provider)) < 2)
    covars_include <- covars_include[-which(covars_include == "data_provider")]
  
  mod_formula <- paste0("value ~ dosage + ",
                        paste0(covars_include, collapse = " + "))
  
  modeled_dat <- polr(formula(mod_formula), data = dat, Hess = T)
  print_res <- tidy(modeled_dat) %>%
    filter(term == "dosage") %>%
    rename(beta = estimate, se = std.error, tstat = statistic) %>%
    dplyr::select(all_of(c("beta", "se", "tstat"))) %>%
    mutate(sample_size = nrow(dat))
  return (print_res)
}

# Apply testing and plotting ----

test_res <- lapply(VARIDS, function (varid) {
  cat(paste0("Running SNP: ", varid, "\n"))
  # Loop through the various adiposity traits
  
  quant_res_tables <- lapply(QUANT_PHENOS, function (qp) {
    # Create full data with dosage and genotype info
    model_dat <- dat_for_replication[[qp]]
    model_dat <- addVarDosage(model_dat, varid)
    
    # Run regressions and plots
    sub_dat <- model_dat %>% filter(visit == 2)
    quant_res <- quantTest(sub_dat, "sex_comb") %>%
      mutate(visit_compared = 2,
             sex_strata = "sex_comb",
             pheno_tested = qp)
    return (quant_res)
  })
  quant_res_tables <- bind_rows(quant_res_tables)
  
  # Loop through categorical traits (only one for now)
  cat_res_tables <- lapply(CAT_PHENOS, function (cp) {
    # Create full data with dosage and genotype info
    model_dat <- dat_for_replication[[cp]]
    model_dat <- addVarDosage(model_dat, varid)
    
    sub_dat <- model_dat %>% filter(visit == 1)
    cat_res <- catTest(sub_dat, "sex_comb") %>%
      mutate(visit_compared = 1,
             sex_strata = "sex_comb",
             pheno_tested = cp)
  })
  cat_res_tables <- bind_rows(cat_res_tables) %>%
    mutate(OR = exp(beta), 
           lci = exp(beta - (1.96*se)), uci = exp(beta + (1.96*se)),
           pval = pt(tstat, df = sample_size))
  
  full_res_tables <- bind_rows(quant_res_tables, cat_res_tables) %>%
    mutate(SNP = varid)
  return (full_res_tables)
})

test_res <- bind_rows(test_res)

write.table(test_res, paste0(resdir, "adipo_change_replication.txt"),
            sep = "\t", quote = F, row.names = F)

