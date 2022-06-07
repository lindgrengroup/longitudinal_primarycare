# Author: Samvida S. Venkatesh
# Date: 12/04/22

library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

# Read files ----

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

mapping_file <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/rename_clusters_map_file.txt",
                           sep = "\t", header = T, stringsAsFactors = F)

cluster_membership <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/",
                      p, "_", sx, "_spline_full_models.txt"),
               sep = "\t", header = T, stringsAsFactors = F)
    
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
}) 
names(cluster_membership) <- PHENOTYPES

# Remap clusters, adding a new column for "mapped_k" ----

new_clusters <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    res <- cluster_membership[[p]][[sx]]
    sub_map <- mapping_file %>% filter(strata == paste0(p, "_", sx))
    
    res$mapped_k <- sub_map$new_k[match(res$assigned_k, 
                                        sub_map$assigned_k)]
    res$eid <- as.character(res$eid)
    
    write.table(res, paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/",
                            p, "_", sx, "_remapped_clusters.txt"),
                sep = "\t", row.names = F, quote = F)
    
    return (res)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
}) 
names(new_clusters) <- PHENOTYPES

# Plot results with new clusters to ensure mapping is correct ----

library(splines)
library(zoo)
library(lme4)

## Read data ----

spline_models <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/spline_models/",
                 p, "_full_model.rds"))
})
names(spline_models) <- PHENOTYPES

raw_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/indiv_qcd_data.rds")[PHENOTYPES]
covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[PHENOTYPES]

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220131_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)

## Wrangle data ----

# Make sure covars have characters for eid
covars <- lapply(PHENOTYPES, function (p) {
  res <- covars[[p]] %>%
    select(all_of(c("eid", "baseline_age", "age_sq"))) %>%
    mutate(eid = as.character(eid))
  return (res)
})
names(covars) <- PHENOTYPES

general_covars <- general_covars %>% 
  select(all_of(c("eid", "sex"))) %>%
  mutate(eid = as.character(eid))

model_dat <- lapply(PHENOTYPES, function (p) {
  # Add in covariates and 't'
  res <- raw_dat[[p]] %>%
    mutate(eid = as.character(eid)) %>%
    left_join(covars[[p]], by = "eid") %>%
    left_join(general_covars, by = "eid") %>%
    mutate(t = age_event - baseline_age)
  
  # Split by sex and add cluster membership
  res_list <- lapply(SEX_STRATA, function (sx) {
    sexsplit <- res
    if (sx != "sex_comb") sexsplit <- sexsplit %>% filter(sex == sx)
    
    sexsplit <- sexsplit %>%
      left_join(new_clusters[[p]][[sx]][, c("eid", "mapped_k")],
                by = "eid") %>%
      filter(!is.na(mapped_k)) %>%
      mutate(cluster = factor(as.character(mapped_k),
                              levels = paste0("k", 1:5)))
    return (sexsplit)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
}) 
names(model_dat) <- PHENOTYPES

## Function to plot observed data on population level, within clusters ----

plotObservedTrajCluster <- function (p, sx) {
  # Round ages to nearest 0.25 yrs
  summ_dat <- model_dat[[p]][[sx]] %>% 
    mutate(age_bin = plyr::round_any(age_event, 0.25, f = round)) %>%
    filter(age_bin >= 30 & age_bin <= 70) %>%
    group_by(cluster, age_bin) %>% 
    summarise(mean_value = mean(value),
              sd_value = sd(value), 
              n = n()) %>%
    mutate(lci_mean = mean_value - 1.96*(sd_value/sqrt(n)),
           uci_mean = mean_value + 1.96*(sd_value/sqrt(n))) 
  
  # Get rolling average mean across 2 years 
  summ_dat <- summ_dat %>% 
    ungroup() %>% 
    group_by(cluster) %>% 
    arrange(age_bin, .by_group = T) %>%
    mutate(interval_width = seq_along(age_bin) - 
             findInterval(age_bin - 2, age_bin),
           mean_value_rolled = rollapply(mean_value, interval_width, mean, 
                                         fill = NA),
           cluster = factor(as.character(cluster),
                            levels = paste0("k", 1:5)))
  
  all_plot <- ggplot(summ_dat,
                     aes(x = age_bin, y = mean_value, 
                         color = cluster, fill = cluster)) +
    # Faintly plot mean and SE in background
    geom_ribbon(aes(ymin = lci_mean, 
                    ymax = uci_mean), alpha = 0.2, 
                linetype = 0) +
    geom_line(aes(y = mean_value), alpha = 0.2) +
    # Add a thicker line for rolling average 
    geom_line(aes(y = mean_value_rolled), alpha = 1) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "Age (years)", 
         y = paste0("Mean and 95% C.I. of mean of ", p), 
         title = paste0("Observed trajectories in each cluster of phenotype: ",
                        p, " strata: ", sx))
  
  return (all_plot)
  
}

## Function to plot fitted data on population level, within clusters ----

plotFittedTrajCluster <- function (p, sx) {
  # Add model fits
  df <- model_dat[[p]][[sx]]
  df$fit <- fitted(spline_models[[p]][[sx]])
  
  # Round ages to nearest 0.25 yrs
  summ_dat <- df %>% 
    mutate(t = plyr::round_any(t, 0.25, f = round)) %>%
    filter(t <= 20) %>%
    group_by(cluster, t) %>% 
    summarise(mean_fit = mean(fit),
              sd_fit = sd(fit), 
              n = n()) %>%
    mutate(lci_mean = mean_fit - 1.96*(sd_fit/sqrt(n)),
           uci_mean = mean_fit + 1.96*(sd_fit/sqrt(n)),
           cluster = factor(as.character(cluster),
                            levels = paste0("k", 1:5))) 
  
  all_plot <- ggplot(summ_dat,
                     aes(x = t, y = mean_fit, 
                         color = cluster, fill = cluster)) +
    geom_ribbon(aes(ymin = lci_mean, 
                    ymax = uci_mean), alpha = 0.2, 
                linetype = 0) +
    geom_line(aes(y = mean_fit)) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "Time from baseline measurement", 
         y = paste0("Mean and 95% C.I. of mean of ", p), 
         title = paste0("Fitted trajectories in each cluster of phenotype: ",
                        p, " strata: ", sx))
  
  return (all_plot)
  
}

## Apply plotting functions ----

plot_list <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    
    pdf(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/plots/trajectories/remapped_",
               p, "_", sx, "_trajectories.pdf"), onefile = T)
    print(plotObservedTrajCluster(p, sx))
    print(plotFittedTrajCluster(p, sx))
    dev.off()
    
    return ()
  })
  return ()
})

# Count overlaps between weight and BMI in newly mapped clusters ----

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

cluster_membership <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/",
                      p, "_", sx, "_remapped_clusters.txt"),
               sep = "\t", header = T, stringsAsFactors = F)
    
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
}) 
names(cluster_membership) <- PHENOTYPES

## Summarise counts for confusion matrix ----

summ_weight_vs_bmi <- lapply(SEX_STRATA, function (sx) {
  wt_ks <- cluster_membership[["Weight"]][[sx]] %>% 
    select(all_of(c("eid", "mapped_k"))) %>%
    rename(weight_k = mapped_k)
  bmi_ks <- cluster_membership[["BMI"]][[sx]] %>% 
    select(all_of(c("eid", "mapped_k"))) %>%
    rename(bmi_k = mapped_k)
  
  to_summ <- inner_join(wt_ks, bmi_ks, by = "eid")
  
  summ_dat <- to_summ %>% group_by(weight_k, bmi_k) %>% 
    count() %>%
    pivot_wider(names_from = bmi_k, values_from = n, values_fill = 0)
  colnames(summ_dat)[-1] <- paste0("bmi_", colnames(summ_dat)[-1])
  
  write.table(summ_dat, 
              paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/logs/weight_bmi_overlap_",
                     sx, ".txt"),
              sep = "\t", row.names = F, quote = F)
  return ()
})

summ_sex_strata <- lapply(PHENOTYPES, function (p) {
  
  F_ks <- cluster_membership[[p]][["F"]] %>% 
    select(all_of(c("eid", "mapped_k"))) %>%
    rename(F_k = mapped_k)
  M_ks <- cluster_membership[[p]][["M"]] %>% 
    select(all_of(c("eid", "mapped_k"))) %>%
    rename(M_k = mapped_k)
  sc_ks <- cluster_membership[[p]][["sex_comb"]] %>% 
    select(all_of(c("eid", "mapped_k"))) %>%
    rename(sc_k = mapped_k)
  
  compare_F <- inner_join(F_ks, sc_ks, by = "eid")
  summ_F <- compare_F %>% group_by(F_k, sc_k) %>% 
    count() %>%
    pivot_wider(names_from = sc_k, values_from = n, values_fill = 0)
  colnames(summ_F)[-1] <- paste0("sex_comb_", colnames(summ_F)[-1])
  write.table(summ_F, 
              paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/logs/",
                     p, "_F_sex_comb_overlap.txt"),
              sep = "\t", row.names = F, quote = F)
  
  compare_M <- inner_join(M_ks, sc_ks, by = "eid")
  summ_M <- compare_M %>% group_by(M_k, sc_k) %>% 
    count() %>%
    pivot_wider(names_from = sc_k, values_from = n, values_fill = 0)
  colnames(summ_M)[-1] <- paste0("sex_comb_", colnames(summ_M)[-1])
  write.table(summ_M, 
              paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/logs/",
                     p, "_M_sex_comb_overlap.txt"),
              sep = "\t", row.names = F, quote = F)
  
  return ()
})

