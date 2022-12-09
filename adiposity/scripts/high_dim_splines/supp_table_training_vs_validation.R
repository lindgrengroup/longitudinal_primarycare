# Author: Samvida S. Venkatesh
# Date: 24/05/2022

library(tidyverse)
theme_set(theme_bw())

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

resdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/clustering/training_vs_validation/"

# Load data ----

# Soft clustering probabilities
soft_clust_res <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    res <- read.table(paste0(resdir, p, "_", sx, "/soft_clustering_probs_", p, "_", sx, ".txt"),
               sep = "\t", header = T, stringsAsFactors = F)
    return (res)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(soft_clust_res) <- PHENOTYPES

# Covariates
covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[PHENOTYPES]
covars <- lapply(covars, function (df) {
  df$eid <- as.character(df$eid)
  return (df)
})

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

# Training vs validation ids
id_class <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    training_ids <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/", 
                                      p, "_", sx, "/ids_training.txt"),
                               sep = "\t", header = F, stringsAsFactors = F)$V1
    training_ids <- data.frame(eid = as.character(training_ids),
                               id_type = "training")
    validation_ids <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/", 
                                        p, "_", sx, "/ids_validation.txt"),
                                 sep = "\t", header = F, stringsAsFactors = F)$V1
    validation_ids <- data.frame(eid = as.character(validation_ids),
                                 id_type = "validation")
    res <- bind_rows(training_ids, validation_ids)
    return (res)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(id_class) <- PHENOTYPES

# Convert soft clustering probabilities to hard assignments ----

hard_clust_res <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    tmp <- soft_clust_res[[p]][[sx]][, -1]
    kclusts <- colnames(tmp)
    assnd_k <- kclusts[apply(tmp, 1, function (x) which.max(x))]
    res <- data.frame(eid = as.character(soft_clust_res[[p]][[sx]]$eid),
                      cluster = assnd_k)
    return (res)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(hard_clust_res) <- PHENOTYPES

# Get summary tables ----

sumtables <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    full_dat <- left_join(hard_clust_res[[p]][[sx]], 
                          id_class[[p]][[sx]], by = "eid")
    full_dat <- left_join(full_dat, covars[[p]], by = "eid")
    full_dat <- left_join(full_dat, general_covars, by = "eid")
    
    sumstats_table <- full_dat %>% 
      group_by(cluster, id_type) %>%
      summarise(count = n(), 
                mean_FU_n = mean(FU_n), sd_FU_n = sd(FU_n),
                to_print_FU_n = paste0(signif(mean_FU_n, 3), " (", signif(sd_FU_n, 3), ")"),
                mean_FUyrs = mean(FUyrs), sd_FUyrs = sd(FUyrs),
                to_print_FUyrs = paste0(signif(mean_FUyrs, 3), " (", signif(sd_FUyrs, 3), ")"),
                mean_bl_age = mean(baseline_age), 
                sd_bl_age = sd(baseline_age),
                to_print_bl_age = paste0(signif(mean_bl_age, 3), " (", signif(sd_bl_age, 3), ")"),
                median_bl_trait = median(baseline_trait), 
                iqr_bl_trait = paste(signif(quantile(baseline_trait, 0.25), 3), 
                                     signif(quantile(baseline_trait, 0.75), 3),
                                     sep = ", "),
                to_print_bl_trait = paste0(signif(median_bl_trait, 3), " (", iqr_bl_trait, ")")) %>%
      ungroup() %>% group_by(id_type) %>%
      mutate(percent = count / sum(count),
             to_print_n = paste0(count, " (", signif(percent*100, 3), "%)")) %>%
      select(all_of(c("cluster", "id_type",
                      "to_print_n", "to_print_FU_n", "to_print_FUyrs", 
                      "to_print_bl_age", "to_print_bl_trait")))
    write.table(sumstats_table,
                paste0(resdir, p, "_", sx, "_cluster_properties.txt"),
                sep = "\t", row.names = F, quote = F)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(sumtables) <- PHENOTYPES
