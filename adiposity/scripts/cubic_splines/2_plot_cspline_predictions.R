# Author: Samvida S. Venkatesh
# Date: 01/11/21

library(lme4)
library(splines)
library(tidyverse)
theme_set(theme_bw())

set.seed(011121)

# Read files ----

PHENOTYPES <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/pheno_names.txt")$V1

slope_models <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/age_adj_",
                 p, ".rds"))
})
names(slope_models) <- PHENOTYPES

blups <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/age_adj_blups_",
                 p, ".rds"))
})
names(blups) <- PHENOTYPES

SEX_STRATA <- c("F", "M")

dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/indiv_qcd_data.rds")[PHENOTYPES]
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/covariates.rds")[PHENOTYPES]

COVARS <- c("baseline_age", "age_sq")

# Wrangle data ----

# Add in covariates
model_dat <- lapply(PHENOTYPES, function (p) {
  res <- merge(dat[[p]], covars[[p]], by = "eid")
  res <- res %>% mutate(t = age_event - baseline_age)
  return (res)
})
names(model_dat) <- PHENOTYPES

# Plot model predictions for selected individuals ----

## Random IDs ----

get_rand_eids <- function (p, sx, n = 25) {
  sub_dat <- covars[[p]] %>% filter(sex == sx)
  nsample <- min(nrow(sub_dat), n)
  sample_eids <- sample(sub_dat$eid, nsample, replace = F)
  return (sample_eids)
}

## IDs with only 2 repeat measures or several repeat measures ----

get_2fu_ids <- function (p, sx, n = 25) {
  sub_dat <- covars[[p]] %>% filter(sex == sx & FU_n == 2)
  nsample <- min(nrow(sub_dat), n)
  sample_eids <- sample(sub_dat$eid, nsample, replace = F)
  return (sample_eids)
}

get_many_fu_ids <- function (p, sx, n = 25) {
  sub_dat <- covars[[p]] %>% filter(sex == sx)
  # Sample from top 10% of FU
  threshold <- quantile(sub_dat$FU_n, 0.9)
  relevant_eids <- sub_dat$eid[sub_dat$FU_n > threshold]
  nsample <- min(length(relevant_eids), n)
  sample_eids <- sample(relevant_eids, nsample, replace = F)
  return (sample_eids)
}

## IDs with very high or very low baseline ages ----

get_tail_baseline_age_ids <- function (p, sx, high = T, n = 25) {
  sub_dat <- covars[[p]] %>% filter(sex == sx)
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
  sub_dat <- covars[[p]] %>% filter(sex == sx)
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

## IDs with very high or very low BLUPs ----

get_tail_BLUP_ids <- function (p, sx, term, high = T, n = 25) {
  # Sample from top or bottom 10% of BLUP
  relevant_dat <- data.frame(eid = rownames(blups[[p]][[sx]]),
                             term_blup = blups[[p]][[sx]][, term])
  if (high) { 
    threshold <- quantile(relevant_dat$term_blup, 0.9) 
    relevant_eids <- relevant_dat$eid[relevant_dat$term_blup > threshold]
  } else { 
    threshold <- quantile(relevant_dat$term_blup, 0.1) 
    relevant_eids <- relevant_dat$eid[relevant_dat$term_blup < threshold]
  }
  nsample <- min(length(relevant_eids), n)
  sample_eids <- sample(relevant_eids, nsample, replace = F)
  return (sample_eids)
}

## Function to create predicted data for set of ids ----

create_prediction_df <- function (p, sx, ids) {
  sub_dat <- subset(model_dat[[p]], 
                    model_dat[[p]]$eid %in% ids)
  # Calculate maximum time-point to predict to for each eid
  sub_dat <- sub_dat %>% group_by(eid) %>% 
    mutate(max_t = max(t))
  
  # Create new data to predict from
  new_data <- sub_dat %>% select(all_of(c("eid", COVARS,
                                          "max_t", "sex"))) %>% 
    distinct(eid, data_provider, .keep_all = T) 
  # Timepoints to extend to 
  t <- unlist(lapply(new_data$max_t, function (x) { 
    seq(0, x, length.out = 30) }))
  tmp <- data.frame(eid = rep(new_data$eid, each = 30),
                    t = t)
  new_data <- merge(tmp, new_data, by = "eid")
  
  # Predict new values
  fitted_results <- 
    as.data.frame(predict(slope_models[[p]][[sx]],
                          newdata = new_data))
  colnames(fitted_results) <- "fit"
  pred_df <- bind_cols(new_data, fitted_results) %>%
    mutate(age_event = baseline_age + t)
  
  return (pred_df)
}

## Function to create plots ----

sex_col_palette <- c("#F8766D", "#00BFC4", "#C77CFF")
names(sex_col_palette) <- c("F", "M", "sex_comb")

plot_predictions <- function (p, sx, ids) {
  raw_dat <- model_dat[[p]] %>% filter(eid %in% ids)
  
  sex_specific <- create_prediction_df(p, sx, ids) %>%
    mutate(model_strata = sx)
  sex_combined <- create_prediction_df(p, "sex_comb", ids) %>%
    mutate(model_strata = "sex_comb")
  
  plot_dat <- bind_rows(sex_specific, sex_combined)
  
  res <- ggplot(plot_dat, aes(x = age_event)) +
    facet_wrap(~eid, ncol = 5) +
    geom_point(data = raw_dat, aes(y = value),
               colour = "black") +
    geom_line(aes(y = fit, colour = model_strata)) +
    scale_color_manual(values = sex_col_palette) +
    labs(x = "Age (years)",
         y = p)
  return (res)
}

## Apply plotting functions ----

plot_list <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(c("F", "M"), function (sx) {
    
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
    
    # Low or high BLUPs for each term
    term_names <- colnames(blups[[p]][[sx]])
    
    low_blup_plots <- lapply(term_names, function (tm_name) {
      res_plot <- plot_predictions(p, sx, 
                                   get_tail_BLUP_ids(p, sx, tm_name, 
                                                     high = F, 25))
      res_plot <- res_plot +
        labs(title = paste0(sx, " in bottom 10% of BLUP for term ", tm_name,
                            " phenotype: ", p))
      return (res_plot)
    })
    
    high_blup_plots <- lapply(term_names, function (tm_name) {
      res_plot <- plot_predictions(p, sx, 
                                   get_tail_BLUP_ids(p, sx, tm_name, 
                                                     high = T, 25))
      res_plot <- res_plot +
        labs(title = paste0(sx, " in top 10% of BLUP for term ", tm_name,
                            " phenotype: ", p))
      return (res_plot)
    })
    
    pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/cubic_splines/predictions/age_adj_",
               p, "_", sx, ".pdf"))
    print(rand_plots)
    print(fu2_plots)
    print(fu_many_plots)
    print(low_bl_age_plots)
    print(high_bl_age_plots)
    print(low_bl_trait_plots)
    print(high_bl_trait_plots)
    print(low_blup_plots)
    print(high_blup_plots)
    dev.off()
    
    return ()
  })
  return ()
})

