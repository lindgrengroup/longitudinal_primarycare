# Author: Samvida S. Venkatesh
# Date: 11/10/21

library(tidyverse)
library(lme4)
library(broom)
library(foreign)
library(MASS)
library(ggpubr)
theme_set(theme_bw())

# Read data ----

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

MAIN_COVARS <- c("baseline_age", "age_sq", "FUyrs", "FU_n",
                 "year_of_birth") # add sex and UKB assessment centre if needed
GEN_COVARS <- paste0("PC", 1:21) # add genotyping array if needed

resdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/"

for_linreg_dat <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SEX_STRATA, function (sx) {
    df <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/replication_genotypes/traits_for_GWAS/",
                            p, "_", sx, ".txt"), 
                     sep = "\t", header = T, stringsAsFactors = F, 
                     comment.char = "$")
    colnames(df)[1] <- "eid"
    df$eid <- as.character(df$eid)
    return (df)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(for_linreg_dat) <- PHENOTYPES

selfrep_wtchg <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/data/main_data_adipo_change.rds")[["Weight_change_1yr"]]

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)
general_covars <- general_covars %>%
  dplyr::select(any_of(c("eid", "sex", "UKB_assmt_centre",
                         MAIN_COVARS, GEN_COVARS))) %>%
  mutate(sex = factor(sex),
         UKB_assmt_centre = factor(UKB_assmt_centre))

apoe_dosage <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sample_variant_counts/rs429358_dosages.txt",
                          sep = " ", header = T, stringsAsFactors = F)
# Remove first row, which contains info on type of column and columns 
# 2, 3, 4 (ID repeat, missingness, sex)
apoe_dosage <- apoe_dosage[-1, c(1, 5)]
colnames(apoe_dosage) <- c("eid", "dosage")

# IDs with AD or dementia
ids_with_dementia <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eids_with_dementia.txt",
                                sep = "\t", header = F, stringsAsFactors = F)$V1

# IDs passed QC for weight-change
selfrep_wtchg_ids_passed_qc <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/adipo_change/sample_qc/Weight_change_1yr_ids_passed_qc.txt",
                                          sep = "\t", header = T, stringsAsFactors = F)
selfrep_wtchg_ids_passed_qc <- selfrep_wtchg_ids_passed_qc %>%
  mutate(eid = as.character(eid),
         genotyping.array = factor(genotyping.array))

# Wrangle data ----

# Convert rs429358 dosage to genotype
DOSAGE_THRESHOLD_0 <- 0.5
DOSAGE_THRESHOLD_2 <- 1.5
apoe_genotypes <- apoe_dosage %>% 
  mutate(genotype = ifelse(dosage < DOSAGE_THRESHOLD_0, "0",
                           ifelse(dosage > DOSAGE_THRESHOLD_2, 
                                  "2", "1")),
         eid = as.character(eid)) %>%
  dplyr::select(all_of(c("eid", "genotype")))

# Add in sex to self-reported weight change and split by sex
# Add in model covariates
# Tally visit number and changes from first visit needed for model
tmp <- selfrep_wtchg %>% 
  group_by(eid) %>% 
  mutate(value = factor(value, 
                        levels = c("Loss", "No change", "Gain"),
                        ordered = T),
         visit = seq(n()),
         age_event_sq = age_event^2) %>%
  left_join(general_covars) %>%
  inner_join(selfrep_wtchg_ids_passed_qc)

selfrep_wtchg <- lapply(SEX_STRATA, function (sx) {
  df <- tmp
  if (sx != "sex_comb") {
    df <- df %>% filter(sex == sx)
  }
  return (df)
})
names(selfrep_wtchg) <- SEX_STRATA

# Function to add hard-called genotype
addGenoGroup <- function (dat) {
  res <- dat
  res$genotype <- apoe_genotypes$genotype[match(res$eid, apoe_genotypes$eid)]
  res$genotype <- as.numeric(res$genotype)
  res <- res[!is.na(res$genotype), ]
  return (res)
}

# Function to add dementia diagnosis
addDementia <- function (dat) {
  res <- dat
  res$dementia <- res$eid %in% ids_with_dementia
  return (res)
}

# Testing functions ----

b1Test <- function (dat, ss) {
  to_ret <- tryCatch({
    covars_include <- GEN_COVARS
    if (ss == "sex_comb") covars_include <- c(covars_include, "sex")
    if (length(unique(dat$genotyping.array)) > 1) covars_include <- c(covars_include, "genotyping.array")
    
    model_formula <- formula(paste0("b1 ~ genotype + ", 
                                    paste0(covars_include, collapse = " + ")))
    
    modeled_dat <- lm(model_formula, data = dat)
    
    print_res <- tidy(modeled_dat) %>%
      filter(term == "genotype") %>%
      rename(beta = estimate, se = std.error, tstat = statistic, pval = p.value) %>%
      dplyr::select(all_of(c("beta", "se", "tstat", "pval"))) %>%
      mutate(sample_size = nrow(dat))
    return (print_res)
  }, error = function (e) {
    return (NULL)
  })
  return (to_ret)
}

clustTest <- function (dat, clust_test, ss) {
  to_ret <- tryCatch({
    covars_include <- c("baseline_trait", MAIN_COVARS, GEN_COVARS)
    if (ss == "sex_comb") covars_include <- c(covars_include, "sex")
    if (length(unique(dat$UKB_assmt_centre)) > 1) covars_include <- c(covars_include, "UKB_assmt_centre")
    if (length(unique(dat$genotyping.array)) > 1) covars_include <- c(covars_include, "genotyping.array")
    
    model_formula <- formula(paste0(clust_test, " ~ genotype + ", 
                                    paste0(covars_include, collapse = " + ")))
    
    modeled_dat <- lm(model_formula, data = dat)
    
    print_res <- tidy(modeled_dat) %>%
      filter(term == "genotype") %>%
      rename(beta = estimate, se = std.error, tstat = statistic, pval = p.value) %>%
      mutate(or = exp(beta), lci = exp(beta-1.96*se), uci = exp(beta+1.96*se)) %>%
      dplyr::select(all_of(c("beta", "se", "or", "lci", "uci",
                             "tstat", "pval"))) %>%
      mutate(sample_size = nrow(dat))
    
    return (print_res)
  }, error = function (e) {
    return (NULL)
  })
  return (to_ret)
}

catTest <- function (dat, ss) {
  to_ret <- tryCatch({
    # Get the correct covariates for adjustment
    # Slightly different set than the quant traits because we shouldn't be
    # using the correlated age-event-sq
    covars_include <- c("BMI", "age_event", "year_of_birth", GEN_COVARS)
    if (ss == "sex_comb") covars_include <- c(covars_include, "sex")
    if (length(unique(dat$UKB_assmt_centre)) > 1) covars_include <- c(covars_include, "UKB_assmt_centre")
    if (length(unique(dat$genotyping.array)) > 1) covars_include <- c(covars_include, "genotyping.array")
    
    mod_formula <- paste0("value ~ genotype + ",
                          paste0(covars_include, collapse = " + "))
    
    modeled_dat <- polr(formula(mod_formula), data = dat, Hess = T)
    
    print_res <- tidy(modeled_dat) %>%
      filter(term == "genotype") %>%
      rename(beta = estimate, se = std.error, tstat = statistic) %>%
      mutate(or = exp(beta), lci = exp(beta-1.96*se), uci = exp(beta+1.96*se),
             sample_size = nrow(dat),
             pval = pt(tstat, df = sample_size)) %>%
      dplyr::select(all_of(c("beta", "se", "or", "lci", "uci",
                             "tstat", "pval", "sample_size"))) 
    return (print_res)
  }, error = function (e) {
    return (NULL)
  })
  return (to_ret)
}

# Apply tests ----

all_res <- lapply(SEX_STRATA, function (sx) {
  cat(paste0("\t", "Sex strata: ", sx, "\n"))
  
  bmi_wt_res <- lapply(PHENOTYPES, function (p) {
    cat(paste0("\t", "b1 phenotype: ", p, "\n"))
    
    full_dat <- addGenoGroup(dat = for_linreg_dat[[p]][[sx]])
    full_dat <- addDementia(dat = full_dat)
    full_dat <- full_dat[complete.cases(full_dat), ]
    
    lmm_full <- b1Test(dat = full_dat, ss = sx) %>%
      mutate(status = "all")
    
    lmm_nd <- b1Test(dat = full_dat[which(!full_dat$dementia), ],
                     ss = sx) %>%
      mutate(status = "no dementia")
    
    lmm_res <- bind_rows(lmm_full, lmm_nd) %>%
      mutate(term = "b1")
    
    clust_res <- lapply(c("k1", "k1_k2", "k1_k2_k3"), function (clustk) {
      cat(paste0("\t", clustk, " phenotype: ", p, "\n"))
      
      clustk_full <- clustTest(dat = full_dat, 
                              clust_test = clustk,
                              ss = sx) %>%
        mutate(status = "all")
      
      clustk_nd <- clustTest(dat = full_dat[which(!full_dat$dementia), ],
                               clust_test = clustk,
                               ss = sx) %>%
        mutate(status = "no dementia")
      
      clustk_res <- bind_rows(clustk_full, clustk_nd) %>%
        mutate(term = clustk)
      return (clustk_res)
    })
    clust_res <- bind_rows(clust_res)
    
    res <- bind_rows(clust_res, lmm_res)
    res$pheno_tested <- p
    return (res)
  })
  bmi_wt_res <- bind_rows(bmi_wt_res)
  
  selfrep_res <- lapply(c(1:3), function (vc) {
    cat(paste0("\t", "Visit: ", vc, " Weight_change_1yr", "\n"))
    sub_dat <- selfrep_wtchg[[sx]] %>% 
      filter(visit == vc)
    sub_dat <- addGenoGroup(dat = sub_dat) 
    sub_dat <- addDementia(dat = sub_dat)
    sub_dat <- sub_dat[complete.cases(sub_dat), ]
    
    res_full <- catTest(dat = sub_dat, ss = sx) %>%
      mutate(status = "all")
    res_nd <- catTest(dat = sub_dat[which(!sub_dat$dementia), ], ss = sx) %>%
      mutate(status = "no dementia")
    
    res <- bind_rows(res_full, res_nd)
    res$visit_compares <- vc
    return (res)
  })
  selfrep_res <- bind_rows(selfrep_res) 
  selfrep_res$pheno_tested <- "Weight_change_1yr"
  
  res <- bind_rows(bmi_wt_res, selfrep_res)
  res$sex_strata <- sx
  return (res)
})
all_res <- bind_rows(all_res)

write.table(all_res, 
            paste0(resdir, "rs429358_replication_by_dementia_status.txt"),
            sep = "\t", quote = F, row.names = F)
