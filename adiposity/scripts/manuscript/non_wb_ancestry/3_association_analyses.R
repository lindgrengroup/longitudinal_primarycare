# Author: Samvida S. Venkatesh
# Date: 31/08/22

library(tidyverse)
library(lme4)
library(broom)
library(foreign)
library(MASS)
library(ggpubr)
theme_set(theme_bw())

# Read in arguments ----

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")
ANCESTRIES <- c("asian", "black", "chinese", "mixed", "other", "white")

MAIN_COVARS <- c("baseline_age", "age_sq", "FUyrs", "FU_n",
                 "year_of_birth") # add sex and UKB assessment centre if needed
GEN_COVARS <- paste0("PC", 1:21) # add genotyping array if needed

infile_path <- "" # REDACTED
gpdat_path <- "" # REDACTED
gen_resources_path <- "" # REDACTED
resdir <- "" # REDACTED

# Read data ----

selfrep_wtchg <- read.table(paste0(infile_path, "/selfrep_wtchg_non_wb.txt"),
                            sep = "\t", header = T, stringsAsFactors = F)
selfrep_wtchg$eid <- as.character(selfrep_wtchg$eid)

blups <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(ANCESTRIES, function (anc) {
    res <- lapply(SEX_STRATA, function (sx) {
      df <- read.table(paste0(infile_path, "/lmm_models/",
                              p, "_", anc, "_", sx, "_all_blups.txt"), 
                       sep = "\t", header = T, stringsAsFactors = F)
      df$eid <- as.character(df$eid)
      return (df)
    })
    names(res) <- SEX_STRATA
    return (res)
  })
  names(res_list) <- ANCESTRIES
  return (res_list)
})
names(blups) <- PHENOTYPES

clustprobs <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(ANCESTRIES, function (anc) {
    res <- lapply(SEX_STRATA, function (sx) {
      df <- read.table(paste0(infile_path, "/highdim_splines/softclust_probs/soft_clustering_probs_",
                              p, "_", anc, "_", sx, ".txt"), 
                       sep = "\t", header = T, stringsAsFactors = F)
      df <- df %>% mutate(eid = as.character(eid),
                          k1_k2 = k1 + k2,
                          k1_k2_k3 = k1 + k2 + k3)
      return (df)
    })
    names(res) <- SEX_STRATA
    return (res)
  })
  names(res_list) <- ANCESTRIES
  return (res_list)
})
names(clustprobs) <- PHENOTYPES

vars_to_replicate <- read.table(paste0(gpdat_path, "/data/lead_snps_to_replicate.txt"),
                                sep = "\t", header = T, stringsAsFactors = F)
VARIDS <- vars_to_replicate$SNP
# For variants without rsids, build variant name
rename_vars <- grep("^chr", VARIDS)
VARIDS[rename_vars] <- paste0(vars_to_replicate$SNP[rename_vars], "_",
                              vars_to_replicate$Tested_Allele[rename_vars], "_", 
                              vars_to_replicate$Other_Allele[rename_vars])
VARIDS <- gsub("^chr", "", VARIDS)

# Genotypes / dosages at variants of interest
var_dosages <- lapply(VARIDS, function (varid) {
  res <- read.table(paste0(gpdat_path, "/sample_variant_counts/",
                           varid, "_dosages.txt"),
                    sep = " ", header = T, stringsAsFactors = F)
  # Remove first row, which contains info on type of column and columns 
  # 2, 3, 4 (ID repeat, missingness, sex)
  res <- res[-1, c(1, 5)]
  colnames(res) <- c("eid", varid)
  return (res)
})
names(var_dosages) <- VARIDS

# Raw data (original) to calculate covariates
raw_dat <- readRDS(paste0(gpdat_path, "/data/non_wb_gp_main_data_passed_longit_filter.rds"))[PHENOTYPES]

general_covars <- read.table(paste0(gen_resources_path, "/220504_QCd_demographic_covariates.txt"),
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)
general_covars <- general_covars %>%
  dplyr::select(any_of(c("eid", "sex", "ancestry", "UKB_assmt_centre",
                         MAIN_COVARS, GEN_COVARS))) %>%
  mutate(sex = factor(sex),
         UKB_assmt_centre = factor(UKB_assmt_centre))

# IDs that passed sample QC
ids_passed_qc <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(ANCESTRIES, function (anc) {
    res <- lapply(SEX_STRATA, function (sx) {
      df <- read.table(paste0(infile_path, "/sample_qc/", 
                              p, "_", anc, "_", sx, "_ids_passed_qc.txt"),
                       sep = "\t", header = T)
      df$eid <- as.character(df$eid)
      return (df)
    })
    names(res) <- SEX_STRATA
    return (res)
  })
  names(res_list) <- ANCESTRIES
  return (res_list)
})
names(ids_passed_qc) <- PHENOTYPES

selfrep_ids_passed_qc <- lapply(ANCESTRIES, function (anc) {
  res <- lapply(SEX_STRATA, function (sx) {
    df <- read.table(paste0(infile_path, "/sample_qc/selfrep_wtchg_", 
                            anc, "_", sx, "_ids_passed_qc.txt"),
                     sep = "\t", header = T)
    df$eid <- as.character(df$eid)
    return (df)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(selfrep_ids_passed_qc) <- ANCESTRIES

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

# Add in sex and ancestry to self-reported weight change and split by these
# Add in model covariates
# Tally visit number and changes from first visit needed for model
tmp <- selfrep_wtchg %>% 
  group_by(eid) %>% 
  mutate(selfrep_wtchg = factor(selfrep_wtchg, 
                                levels = c("Loss", "No change", "Gain"),
                                ordered = T),
         visit = seq(n()),
         age_event_sq = age_event^2) %>%
  left_join(general_covars)

selfrep_wtchg <- lapply(ANCESTRIES, function (anc) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    df <- tmp %>% filter(ancestry == anc)
    if (sx != "sex_comb") {
      df <- df %>% filter(sex == sx)
    }
    return (df)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(selfrep_wtchg) <- ANCESTRIES

covars <- lapply(PHENOTYPES, function (p) {
  df <- raw_dat[[p]] %>%
    group_by(eid) %>% 
    arrange(age_event, .by_group = T) %>%
    summarise(baseline_trait = first(value),
              baseline_age = first(age_event),
              age_sq = baseline_age^2,
              FU_n = n(),
              FUyrs = last(age_event) - first(age_event))
  df <- left_join(df, 
                  general_covars,
                  by = "eid")
  return (df)
})
names(covars) <- PHENOTYPES

# Wrangling functions ----

# Function to add hard-called genotype
addGenoGroup <- function (dat, varid) {
  res <- dat
  var_dat <- var_dosages_hardcall[[varid]]
  res$genotype <- var_dat$genotype[match(res$eid, var_dat$eid)]
  res$genotype <- as.numeric(res$genotype)
  res <- res[!is.na(res$genotype), ]
  return (res)
}

squeezeProbs <- function (x, nsamples = 100) {
  squeezed_x <- (x*(nsamples - 1) + 0.5)/nsamples
  return (squeezed_x)
}

getLogit <- function (x) {
  return (log(x/(1-x)))
}

# Testing functions ----

getLMMTerm <- function (p, anc, ss) {
  lmm_df <- blups[[p]][[anc]][[ss]]
  # Only retain the IDs that passed genotyping sample QC
  lmm_df <- lmm_df %>% filter(eid %in% ids_passed_qc[[p]][[anc]][[ss]]$eid)
  # Add covariates
  lmm_df <- left_join(lmm_df, covars[[p]], by = "eid") 
  
  # Linearly adjust slope for covariates
  covars_include <- MAIN_COVARS
  if (ss == "sex_comb") covars_include <- c(covars_include, "sex")
  if (length(unique(lmm_df$UKB_assmt_centre)) > 1) covars_include <- c(covars_include, "UKB_assmt_centre")
  
  model_formula <- formula(paste0("t ~ X.Intercept. + ", 
                                  paste0(covars_include, collapse = " + ")))
  
  modelled_slope <- lm(model_formula, lmm_df)
  # RINT
  trait_to_rint <- resid(modelled_slope)
  rinted_trait <- qnorm((rank(trait_to_rint) - 0.5) / sum(!is.na(trait_to_rint)))
  
  # Lmm phenotypes
  res <- data.frame(eid = lmm_df$eid,
                    b1 = rinted_trait)
  return (res)
}

getClustTerm <- function (p, anc, ss) {
  logit_probs <- as.data.frame(apply(clustprobs[[p]][[anc]][[ss]][, c("k1", "k1_k2", "k1_k2_k3")], 2, 
                                     FUN = function (x) {
                                       getLogit(squeezeProbs(x, nsamples = 100))
                                     }))
  logit_probs$eid <- clustprobs[[p]][[anc]][[ss]]$eid
  return (logit_probs)
}

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
    
    mod_formula <- paste0("selfrep_wtchg ~ genotype + ",
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

all_res <- lapply(VARIDS, function (varid) {
  cat(paste0("Running SNP: ", varid, "\n"))
  res_list <- lapply(ANCESTRIES, function (anc) {
    per_sex <- lapply(SEX_STRATA, function (sx) {
      cat(paste0("\t", "Ancestry and sex strata: ", anc, "_", sx, "\n"))
      bmi_wt_res <- lapply(PHENOTYPES, function (p) {
        cat(paste0("\t", "b1 phenotype: ", p, "\n"))
        for_lmm <- addGenoGroup(dat = getLMMTerm(p = p, anc = anc, ss = sx), 
                                varid = varid)
        for_lmm <- for_lmm %>% 
          left_join(covars[[p]]) %>%
          left_join(ids_passed_qc[[p]][[anc]][[sx]])
        
        lmm_res <- b1Test(dat = for_lmm[complete.cases(for_lmm), ], 
                          ss = sx)
        if (!is.null(lmm_res)) lmm_res$term <- "b1"
        
        clust_res <- lapply(c("k1", "k1_k2", "k1_k2_k3"), function (clustk) {
          cat(paste0("\t", clustk, " phenotype: ", p, "\n"))
          for_clustk <- addGenoGroup(dat = getClustTerm(p = p, anc = anc, ss = sx),
                                     varid = varid)
          for_clustk <- for_clustk %>% 
            left_join(covars[[p]]) %>%
            left_join(ids_passed_qc[[p]][[anc]][[sx]])
          clustk_res <- clustTest(dat = for_clustk[complete.cases(for_clustk), ], 
                                  clust_test = clustk,
                                  ss = sx) %>%
            mutate(term = clustk)
          return (clustk_res)
        })
        clust_res <- bind_rows(clust_res)
        
        res <- bind_rows(clust_res, lmm_res)
        if (!is.null(res)) res$pheno_tested <- p
        return (res)
      })
      bmi_wt_res <- bind_rows(bmi_wt_res)
      
      selfrep_res <- lapply(c(1:3), function (vc) {
        cat(paste0("\t", "Visit: ", vc, " Weight_change_1yr", "\n"))
        sub_dat <- selfrep_wtchg[[anc]][[sx]] %>% 
          filter(visit == vc)
        sub_dat <- addGenoGroup(dat = sub_dat, 
                                varid = varid) %>%
          left_join(selfrep_ids_passed_qc[[anc]][[sx]])
        res <- catTest(dat = sub_dat[complete.cases(sub_dat), ], 
                       ss = sx) 
        if(!is.null(res)) res$visit_compared <- vc
        return (res)
      })
      selfrep_res <- bind_rows(selfrep_res) 
      if (!is.null(selfrep_res)) selfrep_res$pheno_tested <- "Weight_change_1yr"
      
      res <- bind_rows(bmi_wt_res, selfrep_res)
      res$sex_strata <- sx
      return (res)
    })
    per_sex <- bind_rows(per_sex)
    per_sex$ancestry <- anc
    return (per_sex)
  })
  res_list <- bind_rows(res_list)
  res_list$SNP <- varid
  return (res_list)
})
all_res <- bind_rows(all_res)

write.table(all_res, 
            paste0(resdir, "all_ancestries_all_phenos.txt"),
            sep = "\t", quote = F, row.names = F)
