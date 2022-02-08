# Author: Samvida S. Venkatesh
# Date: 01/11/21

library(lme4)
library(splines)
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

set.seed(011121)

# Read files ----

PHENOTYPES <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/pheno_names.txt")$V1

lmm_models <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/lmm/lmm_",
                 p, ".rds"))
})
names(lmm_models) <- PHENOTYPES

cspline_models <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/age_adj_",
                 p, ".rds"))
})
names(cspline_models) <- PHENOTYPES

SEX_STRATA <- c("F", "M", "sex_comb")

dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/indiv_qcd_data.rds")[PHENOTYPES]
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/covariates.rds")[PHENOTYPES]

PCs <- paste0("PC", 1:21)
COVARS <- c("baseline_age", "age_sq", "FUyrs", "data_provider")

# Wrangle data ----

# Get time from baseline measurement and add in covariates
model_dat <- lapply(PHENOTYPES, function (p) {
  res <- merge(dat[[p]], covars[[p]], by = "eid")
  res <- res %>% mutate(t = age_event - baseline_age)
  return (res)
})
names(model_dat) <- PHENOTYPES

# Plot model predictions for selected individuals ----

## Random IDs ----

get_rand_eids <- function (p, sx, n = 25) {
  sub_dat <- covars[[p]]
  if (sx != "sex_comb") {
    sub_dat <- sub_dat %>% filter(sex == sx)
  }
  nsample <- min(nrow(sub_dat), n)
  sample_eids <- sample(sub_dat$eid, nsample, replace = F)
  return (sample_eids)
}

## IDs with only 2 repeat measures or several repeat measures ----

get_2fu_ids <- function (p, sx, n = 25) {
  sub_dat <- covars[[p]] %>% filter(FU_n == 2)
  if (sx != "sex_comb") {
    sub_dat <- sub_dat %>% filter(sex == sx)
  }
  nsample <- min(nrow(sub_dat), n)
  sample_eids <- sample(sub_dat$eid, nsample, replace = F)
  return (sample_eids)
}

get_many_fu_ids <- function (p, sx, n = 25) {
  sub_dat <- covars[[p]]
  if (sx != "sex_comb") {
    sub_dat <- sub_dat %>% filter(sex == sx)
  }
  # Sample from top 10% of FU
  threshold <- quantile(sub_dat$FU_n, 0.9)
  relevant_eids <- sub_dat$eid[sub_dat$FU_n > threshold]
  nsample <- min(length(relevant_eids), n)
  sample_eids <- sample(relevant_eids, nsample, replace = F)
  return (sample_eids)
}

## IDs with very high or very low baseline ages ----

get_tail_baseline_age_ids <- function (p, sx, high = T, n = 25) {
  sub_dat <- covars[[p]]
  if (sx != "sex_comb") {
    sub_dat <- sub_dat %>% filter(sex == sx)
  }
  # Sample from top or bottom 10% of baseline age
  if (high) { 
    threshold <- quantile(sub_dat$baseline_age, 0.9) 
    relevant_eids <- sub_dat$eid[sub_dat$baseline_age > threshold]
  } else { 
    threshold <- quantile(sub_dat$baseline_age, 0.1)
    relevant_eids <- sub_dat$eid[sub_dat$baseline_age < threshold]
  }
  nsample <- min(length(relevant_eids), n)
  sample_eids <- sample(relevant_eids, nsample, replace = F)
  return (sample_eids)
}

## IDs with very high or very low baseline trait ----

get_tail_baseline_trait_ids <- function (p, sx, high = T, n = 25) {
  sub_dat <- covars[[p]]
  if (sx != "sex_comb") {
    sub_dat <- sub_dat %>% filter(sex == sx)
  }
  # Sample from top or bottom 10% of baseline age
  if (high) { 
    threshold <- quantile(sub_dat$baseline_trait, 0.9) 
    relevant_eids <- sub_dat$eid[sub_dat$baseline_trait > threshold]
  } else { 
    threshold <- quantile(sub_dat$baseline_trait, 0.1)
    relevant_eids <- sub_dat$eid[sub_dat$baseline_trait < threshold]
  }
  nsample <- min(length(relevant_eids), n)
  sample_eids <- sample(relevant_eids, nsample, replace = F)
  return (sample_eids)
}

## Function to create predicted data for set of ids ----

create_pred_df <- function (model_type = "lmm", p, sx, ids) {
  sub_dat <- subset(model_dat[[p]], 
                    model_dat[[p]]$eid %in% ids)
  # Calculate maximum time-point to predict to for each eid
  sub_dat <- sub_dat %>% group_by(eid) %>% 
    mutate(max_t = max(t))
  
  # Create new data to predict from
  new_data <- sub_dat %>% select(all_of(c("eid", COVARS, PCs,
                                          "max_t", "sex"))) %>% 
    distinct(eid, data_provider, .keep_all = T) 
  # Timepoints to extend to 
  t <- unlist(lapply(new_data$max_t, function (x) { 
    seq(0, x, length.out = 30) }))
  tmp <- data.frame(eid = rep(new_data$eid, each = 30),
                    data_provider = rep(new_data$data_provider, each = 30),
                    t = t)
  new_data <- merge(tmp, new_data, by = c("eid", "data_provider"))
  
  # Predict new values
  if (model_type == "lmm") {
    fitted_results <- 
      as.data.frame(predict(lmm_models[[p]][[sx]],
                            newdata = new_data))
  } else if (model_type == "cspline") {
    fitted_results <- 
      as.data.frame(predict(cspline_models[[p]][[sx]],
                            newdata = new_data))
  } 
  
  colnames(fitted_results) <- "fit"
  pred_df <- bind_cols(new_data, fitted_results)
  
  # At each time-point, average across data providers for each individual
  # Add in data on age at predicted event for plot
  pred_df <- pred_df %>% group_by(eid, t) %>%
    summarise(fit = mean(fit)) %>%
    left_join(covars[[p]], by = "eid") %>%
    mutate(age_event = t + baseline_age) %>%
    select(eid, age_event, fit)
  
  return (pred_df)
}

## Function to create plots ----

plot_predictions <- function (p, sx, ids) {
  raw_dat <- model_dat[[p]] %>% filter(eid %in% ids)
  
  lmm_preds <- create_pred_df("lmm", p, sx, ids) %>%
    mutate(model_type = "lmm")
  cspline_preds <- create_pred_df("cspline", p, sx, ids) %>%
    mutate(model_type = "cspline")
  plot_dat <- bind_rows(lmm_preds, cspline_preds)
  
  res <- ggplot(plot_dat, aes(x = age_event)) +
    facet_wrap(~eid, ncol = 5, scales = "free") +
    geom_point(data = raw_dat, aes(y = value),
               colour = "grey", size = 0.5) +
    geom_line(aes(y = fit, colour = model_type)) +
    scale_color_brewer(palette = "Dark2", guide = F) +
    labs(x = "Age (years)",
         y = p)
  
  return (res)
}

## Apply plotting functions ----

plot_list <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    
    # Random ids
    rand_plots <- plot_predictions(p, sx, get_rand_eids(p, sx, 25))
    rand_plots <- rand_plots + 
      labs(title = paste0("Randomly selected ", sx, ", phenotype: ", p))
    
    # Only 2 repeat measures
    fu2_plots <- plot_predictions(p, sx, get_2fu_ids(p, sx, 25))
    fu2_plots <- fu2_plots + 
      labs(title = paste0(sx, " with only 2 measures, phenotype: ", p))
    
    # Many repeat measures
    fu_many_plots <- plot_predictions(p, sx, get_many_fu_ids(p, sx, 25))
    fu_many_plots <- fu_many_plots + 
      labs(title = paste0(sx, " in top 10% of # repeat measures, phenotype: ", 
                          p))
    
    # Low baseline age
    low_bl_age_plots <- 
      plot_predictions(p, sx, get_tail_baseline_age_ids(p, sx, high = F, 25))
    low_bl_age_plots <- low_bl_age_plots + 
      labs(title = paste0(sx, " in bottom 10% of baseline age, phenotype: ", 
                          p))
    
    # High baseline age
    high_bl_age_plots <- 
      plot_predictions(p, sx, get_tail_baseline_age_ids(p, sx, high = T, 25))
    high_bl_age_plots <- high_bl_age_plots + 
      labs(title = paste0(sx, " in top 10% of baseline age, phenotype: ", 
                          p))
    
    # Low baseline trait
    low_bl_trait_plots <- 
      plot_predictions(p, sx, get_tail_baseline_trait_ids(p, sx, high = F, 25))
    low_bl_trait_plots <- low_bl_trait_plots +
      labs(title = paste0(sx, " in bottom 10% of baseline trait, phenotype: ", 
                          p))
    
    # High baseline trait
    high_bl_trait_plots <- 
      plot_predictions(p, sx, get_tail_baseline_trait_ids(p, sx, high = T, 25))
    high_bl_trait_plots <- high_bl_trait_plots + 
      labs(title = paste0(sx, " in top 10% of baseline trait, phenotype: ", 
                          p))
    
    pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/all_model_predictions_",
               p, "_", sx, ".pdf"))
    print(rand_plots)
    print(fu2_plots)
    print(fu_many_plots)
    print(low_bl_age_plots)
    print(high_bl_age_plots)
    print(low_bl_trait_plots)
    print(high_bl_trait_plots)
    dev.off()
    
    return ()
  })
  return ()
})

