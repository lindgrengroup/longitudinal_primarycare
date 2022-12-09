# Author: Samvida S. Venkatesh
# Date: 07/03/2022

# Load an R module with Bioconductor to automatically access these packages
# Otherwise specify your own package libraries to find packages

library(lme4)
library(splines)
library(argparse)
library(tidyverse)
library(RColorBrewer)
library(zoo)
theme_set(theme_bw())

set.seed(070322)

# Parse arguments ----

parser <- ArgumentParser()
parser$add_argument("--idFile", required=TRUE,
                    help = "Two-column whitespace separated .txt file with id and classification; first line is a header")
parser$add_argument("--biomarkers", 
                    default = "all",
                    help = "Biomarkers to plot")
parser$add_argument("--outPrefix", required=TRUE, 
                    help = "Prefix for output logs and plots")
args <- parser$parse_args()
print(args)

# Read data ----

logFile <- paste0(args$outPrefix, ".txt")

# UKB eids and their classification
id_class <- read.table(args$idFile, header = T, stringsAsFactors = F,
                       quote = "", comment.char = "~")
colnames(id_class) <- c("eid", "class")
CLASSES <- sort(unique(id_class$class))

sink(logFile, append = T)
cat(paste0("Number of ids provided: ", nrow(id_class), "\n"))
cat(paste0("Classification: ", "\n"))
print(table(id_class$class))
cat("\n")
sink()

# Biomarkers longitudinal GP data
BIOMS <- strsplit(args$biomarkers, ", |,|\n|\\s")[[1]]
if (BIOMS == "all") {
  BIOMS <- c("BMI", "WC", "Weight", "WHR")
}
biomarker_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/indiv_qcd_data.rds")[BIOMS]
covar_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[BIOMS]

# Cubic spline models (for now only keep sex-combined for simplicity)
models <- lapply(BIOMS, function (bm) {
  res <- readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/age_adj_",
                                 bm, ".rds"))[["sex_comb"]]
  return (res)
})
names(models) <- BIOMS

# Wrangle data ----

# Subset GP data to ids that are in the provided id file and 
# merge with classification
plot_dat <- lapply(BIOMS, function (bm) {
  df <- biomarker_dat[[bm]] %>% ungroup()
  res <- df %>% filter(eid %in% id_class$eid)
  # Merge in covariates for modelling
  res <- left_join(res, covar_dat[[bm]])
  # Merge in classification
  res$class <- id_class$class[match(res$eid, id_class$eid)]
  res <- res %>% mutate(class = factor(as.character(class),
                                       levels = CLASSES),
                        t = age_event - baseline_age,
                        age_event_sq = age_event^2)
  return (res)
})
names(plot_dat) <- BIOMS

# Subset covariate data to ids that are in the provided id file and
# merge with classification
covar_plot_dat <- lapply(BIOMS, function (bm) {
  df <- covar_dat[[bm]] %>% ungroup()
  res <- df %>% filter(eid %in% id_class$eid)
  # Merge in classification
  res$class <- id_class$class[match(res$eid, id_class$eid)]
  res <- res %>% mutate(class = factor(as.character(class),
                                       levels = CLASSES))
  # Subset to relevant data and pivot longer for plot
  res <- res %>% select(eid, class, 
                        baseline_age, baseline_trait,
                        FUyrs, FU_n) %>%
    pivot_longer(cols = all_of(c("baseline_age", "baseline_trait",
                                 "FUyrs", "FU_n")),
                 names_to = "covariate", values_to = "covar_value")
  return (res)
})
names(covar_plot_dat) <- BIOMS

# Count number of individuals in each class to make sampling quicker later
summ_dat <- lapply(BIOMS, function (bm) {
  res <- plot_dat[[bm]] %>% distinct(eid, class) %>%
    group_by(class) %>% mutate(n_group = n())
  
  # Log to file
  sink(logFile, append = T)
  cat(paste0("\t", "** BIOMARKER **", bm, "\n",
             "\t", "Number of ids retained for plot: ", nrow(res), "\n"))
  cat(paste0("\t", "Classification: "))
  print(table(res$class))
  cat("\n")
  sink()
  
  return (res)
})
names(summ_dat) <- BIOMS

# Create color palette based on data classification
if (length(CLASSES) > 8) {
  custom_col_pal <- colorRampPalette(brewer.pal(8, "Set1"))(length(CLASSES))
} else {
  custom_col_pal <- brewer.pal(length(CLASSES), "Set1")
}
names(custom_col_pal) <- CLASSES

# Functions to plot covariate distributions within each class ----

plot_covar_distribution <- function (bm) {
  to_plot <- covar_plot_dat[[bm]] 
  res_plot <- ggplot(to_plot, 
                      aes(x = class, y = covar_value)) +
    facet_wrap(~covariate, ncol = 2, scales = "free_y") +
    geom_violin(aes(fill = class), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_manual(values = custom_col_pal) + 
    labs(x = "Classification", y = "Covariate value", 
         title = paste0("Covariates of ", bm, " distribution in classes")) +
    theme(legend.position = "none")
  return (res_plot)
}

# Functions to create predicted data for ids ----

create_prediction_df <- function (bm, id_df) {
  # Get biomarker data and subset to relevant ids
  sub_dat <- plot_dat[[bm]] %>% filter(eid %in% id_df$eid)
  
  # Calculate maximum time-point to predict to for each eid
  sub_dat <- sub_dat %>% group_by(eid) %>% 
    mutate(max_t = ceiling(max(t)))
  
  # Create new data to predict from
  new_data <- sub_dat %>% 
    select(all_of(c("eid", "baseline_age", "age_sq", "max_t", "sex"))) %>% 
    distinct(eid, .keep_all = T) 
  # Timepoints to extend to
  ts <- lapply(1:nrow(new_data), FUN = function (i) { 
    seq(0, new_data$max_t[i], by = 0.25) })
  length_ts <- unlist(lapply(ts, function (x) length(x)))
  t <- unlist(ts)
  tmp <- data.frame(eid = rep(new_data$eid, length_ts),
                    t = t)
  new_data <- merge(tmp, new_data, by = "eid")
  
  # Predict new values
  fitted_results <- 
    as.data.frame(predict(models[[bm]], newdata = new_data))
  colnames(fitted_results) <- "fit"
  pred_df <- bind_cols(new_data, fitted_results)
  
  # Add column for age at event and class
  pred_df$class <- id_df$class[match(pred_df$eid, id_df$eid)]
  pred_df <- pred_df %>% 
    mutate(age_event = t + baseline_age,
           class = factor(as.character(class),
                          levels = CLASSES)) %>%
    select(eid, class, t, age_event, fit) 
  
  return (pred_df)
}

# Functions to plot random individuals within each class ----

## Get random ids within each class ----

get_rand_eids <- function (bm, n_each = 5) {
  # Return dataframe of id and class
  sampled_ids <- summ_dat[[bm]] %>%
    group_by(class) %>%
    sample_n(ifelse(n_group < n_each, n_group, n_each)) %>%
    select(eid, class)
  return (sampled_ids)
}

## Functions to create plots ----

# One sub-plot per individual, faceted by id
plot_indiv_observed_dat <- function (bm, id_df) {
  # Get biomarker data to plot and subset to relevant ids
  raw_dat <- plot_dat[[bm]] %>% filter(eid %in% id_df$eid)
  
  # Get predicted data to plot
  pred_dat <- create_prediction_df(bm, id_df)
  
  # Order the ids by class for plot
  id_df <- raw_dat %>% distinct(eid, class)
  id_levels <- id_df$eid[order(id_df$class)]
  raw_dat <- raw_dat %>% mutate(eid_f = factor(as.character(eid),
                                               levels = id_levels))
  pred_dat <- pred_dat %>% mutate(eid_f = factor(as.character(eid),
                                               levels = id_levels))
  
  res <- ggplot(pred_dat, aes(x = age_event)) +
    facet_wrap(.~eid_f, ncol = 5, scales = "free_x") +
    geom_point(data = raw_dat, aes(y = value),
               colour = "black") +
    geom_line(aes(y = fit, colour = class)) +
    scale_color_manual(values = custom_col_pal) +
    labs(x = "Age (years)",
         y = bm,
         title = paste0("Randomly selected ids within each class; biomarker: ", bm))
  return (res)
}

# Sample up to 100 individuals in each class and plot them all (faceted by class)
plot_many_observed_dat <- function (bm) {
  # Get up to 100 random IDs within each cluster
  ids_to_plot <- get_rand_eids(bm, n_each = 100)
  # Within each cluster, choose up to 5 random ids to highlight
  highlight_ids <- ids_to_plot %>% group_by(class) %>% 
    mutate(n_group = n()) %>%
    sample_n(ifelse(n_group < 5, n_group, 5))
  highlight_ids <- unique(highlight_ids$eid)
  
  # Get data
  to_plot <- create_prediction_df(bm, ids_to_plot)
  to_plot <- to_plot %>% mutate(highlight = eid %in% highlight_ids)
  
  # Plot
  res_plot <- ggplot(to_plot %>% filter(!highlight), 
                     aes(x = age_event, y = fit, group = eid,
                         colour = class)) +
    facet_wrap(~class) +
    geom_line(alpha = 0.1) +
    geom_line(data = to_plot %>% filter(highlight),
              alpha = 1) +
    scale_color_manual(values = custom_col_pal) +
    labs(x = "Age (years)",
         y = bm,
         title = paste0("Sample of ", bm, " trajectories in each class"))
  return (res_plot)
}

# Functions to plot group mean trajectories within each class ----

## Age binned summaries ----

# Given a dataframe with time-points, add age bins cut by specified interval length,
# summarise values across individuals in each age bin
# and get the rolling average to plot

summ_pred_dat <- function (df, roll_years = 2) {
  res <- df %>% 
    group_by(class, t) %>% 
    summarise(mean_value = mean(fit),
              sd_value = sd(fit), 
              n = n()) %>%
    mutate(lci_mean = mean_value - 1.96*(sd_value/sqrt(n)),
           uci_mean = mean_value + 1.96*(sd_value/sqrt(n))) 
  
  fn_roll <- function (x) { rollapply(x, roll_years/0.25, 
                                      mean, fill = NA) }
  res <- res %>% ungroup() %>% 
    group_by(class) %>% 
    arrange(t, .by_group = T) %>%
    mutate(across(c(mean_value, lci_mean, uci_mean),
                  fn_roll))
  return (res)
}

## Function to create plots ----

plot_mean_obs_dat <- function (bm) {
  # Get data
  # Create predicted df for all ids
  pred_dat <- create_prediction_df(bm, summ_dat[[bm]])
  
  # Apply rolling mean across 2 years
  to_plot <- summ_pred_dat(pred_dat, roll_years = 2) %>%
    filter(t <= 20)
  
  res_plot <- ggplot(to_plot,
                     aes(x = t, y = mean_value,
                         fill = class,
                         colour = class)) +
    geom_line() +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.2) +
    scale_color_manual(values = custom_col_pal) +
    scale_fill_manual(values = custom_col_pal) +
    labs(x = "Time from baseline measurement (years)",
         y = paste0("Mean and 95% C.I. of mean ", bm),
         title = paste0("Modelled trajectories of ", bm, " in each class"))
  return (res_plot)
}

# Run through each biomarker and print plots ----

# Go through each biomarker
lapply(BIOMS, function (bm) {
  pdf(paste0(args$outPrefix, "_", bm, "_trajectories.pdf"), onefile = T)
  # Covariate distribution
  print(plot_covar_distribution(bm))
  # Random sample (individual ids)
  print(plot_indiv_observed_dat(bm, get_rand_eids(bm, 5)))
  # Random sample (100s of ids)
  print(plot_many_observed_dat(bm))
  # Group mean 
  print(plot_mean_obs_dat(bm))
  dev.off()
})
