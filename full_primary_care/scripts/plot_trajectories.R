# Author: Samvida S. Venkatesh
# Date: 07/03/2022

# Load an R module with Bioconductor to automatically access these packages
# Otherwise specify your own package libraries to find packages

library(argparse)
library(tidyverse)
library(RColorBrewer)
library(zoo)
theme_set(theme_bw())

# Parse arguments ----

parser <- ArgumentParser()
parser$add_argument("--idFile", required=TRUE,
                    help = "Two-column whitespace separated .txt file with id and classification; first line is a header")
parser$add_argument("--biomarkers", 
                    default = "all",
                    help = "Biomarkers to plot")
parser$add_argument("--logFile", required=TRUE, 
                    help = "Name of (.txt) file to store logs on number of individuals plotted")
parser$add_argument("--outPlotDir", required=TRUE, 
                    help = "Path to directory to save output plots")
args <- parser$parse_args()
print(args)

# Read data ----

# UKB eids and their classification
id_class <- read.table(args$idFile, header = T, stringsAsFactors = F)
colnames(id_class) <- c("eid", "class")
CLASSES <- sort(unique(id_class$class))

sink(args$logFile, append = T)
cat(paste0("Number of ids provided: ", nrow(id_class), "\n"))
cat(paste0("Classification: ", "\n"))
print(table(id_class$class))
cat("\n")
sink()

# Biomarkers longitudinal GP data
BIOMS <- strsplit(args$biomarkers, ", |,|\n|\\s")[[1]]
if (BIOMS == "all") {
  BIOMS <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/code_lists/qcd_traits_available.txt",
                      header = F, stringsAsFactors = F)$V1
}
biomarker_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/indiv_qcd_data.rds")[BIOMS]
covar_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[BIOMS]

# Wrangle data ----

# Subset GP data to ids that are in the provided id file and 
# merge with classification
plot_dat <- lapply(BIOMS, function (bm) {
  df <- biomarker_dat[[bm]] %>% ungroup()
  res <- df %>% filter(eid %in% id_class$eid)
  # Merge in classification
  res$class <- id_class$class[match(res$eid, id_class$eid)]
  res <- res %>% mutate(class = factor(as.character(class),
                                       levels = CLASSES))
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
  sink(args$logFile, append = T)
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
  to_plot <- plot_dat[[bm]] %>% filter(eid %in% id_df$eid)
  
  # Order the ids by class for plot
  id_df <- to_plot %>% distinct(eid, class)
  id_levels <- id_df$eid[order(id_df$class)]
  to_plot <- to_plot %>% mutate(eid_f = factor(as.character(eid),
                                               levels = id_levels))
  
  res <- ggplot(to_plot, aes(x = age_event, y = value,
                             group = eid)) +
    facet_wrap(.~eid_f, ncol = 5, scales = "free_x") +
    geom_point(aes(y = value), colour = "black") +
    geom_line(aes(colour = class)) +
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
  to_plot <- plot_dat[[bm]] %>% filter(eid %in% ids_to_plot$eid) %>%
    mutate(highlight = eid %in% highlight_ids)
  
  # Plot
  res_plot <- ggplot(to_plot %>% filter(!highlight), 
                     aes(x = age_event, y = value, group = eid,
                         colour = class)) +
    facet_wrap(~class) +
    geom_point(alpha = 0.1) +
    geom_line(alpha = 0.1) +
    geom_point(data = to_plot %>% filter(highlight),
               alpha = 1) +
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

# Given a dataframe with ages, add age bins cut by specified interval length,
# summarise values across individuals in each age bin
# and get the rolling average to plot

summ_by_age_bins <- function (df, interval_yrs = 0.5, roll_years = 5) {
  # Cut ages into 6-month (0.5 yr) bins
  cut_points <- seq(20, 80, by = interval_yrs)
  cut_labels <- seq(20, 80 - interval_yrs, by = interval_yrs)
  df$age_bin <- cut(df$age_event, 
                    breaks = cut_points, 
                    include.lowest = T,
                    labels = cut_labels)
  res <- df %>% 
    group_by(class, age_bin) %>% 
    summarise(mean_value = mean(value),
              sd_value = sd(value), 
              n = n()) %>%
    mutate(lci_mean = mean_value - 1.96*(sd_value/sqrt(n)),
           uci_mean = mean_value + 1.96*(sd_value/sqrt(n))) 
  
  fn_roll <- function (x) { rollapply(x, roll_years/interval_yrs, 
                                      mean, fill = NA) }
  res <- res %>% ungroup() %>% 
    group_by(class) %>% 
    mutate(age_bin = as.numeric(as.character(age_bin))) %>%
    arrange(age_bin, .by_group = T) %>%
    mutate(across(c(mean_value, lci_mean, uci_mean),
                  fn_roll))
  return (res)
}

## Function to create plots ----

plot_mean_obs_dat <- function (bm) {
  # Get data
  to_plot <- summ_by_age_bins(plot_dat[[bm]], interval_yrs = 0.25, roll_years = 5)
  
  res_plot <- ggplot(to_plot,
                     aes(x = age_bin, y = mean_value,
                         fill = class,
                         colour = class)) +
    geom_line() +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.2) +
    scale_color_manual(values = custom_col_pal) +
    scale_fill_manual(values = custom_col_pal) +
    labs(x = "Age (years)",
         y = paste0("Mean and 95% C.I. of mean ", bm),
         title = paste0("Observed trajectories of ", bm, " in each class"))
  return (res_plot)
}

# Run through each biomarker and print plots ----

# Go through each biomarker
lapply(BIOMS, function (bm) {
  pdf(paste0(args$outPlotDir, bm, "_trajectories.pdf"), onefile = T)
  # Covariate distribution
  print(plot_covar_distribution(bm))
  # Random sample (individual ids)
  print(plot_indiv_observed_dat(bm, get_rand_eids(bm, 5)))
  # Random sample (100s of ids)
  print(plot_many_observed_dat(bm))
  # Group mean 
  print(plot_mean_obs_dat(bm))
  # Covariates of interest
  dev.off()
})
