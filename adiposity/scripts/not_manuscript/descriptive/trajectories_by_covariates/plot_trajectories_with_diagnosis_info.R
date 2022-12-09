# Author: Samvida S. Venkatesh
# Date: 07/03/2022

# Load an R module with Bioconductor to automatically access these packages
# Otherwise specify your own package libraries to find packages

library(argparse)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(zoo)
theme_set(theme_bw())

# Parse arguments ----

parser <- ArgumentParser()
parser$add_argument("--diagnoses", 
                    default = "all",
                    help = "Diagnoses to include in the plots")
parser$add_argument("--biomarker", 
                    default = "Weight",
                    help = "Name of biomarker (ONLY ONE) to plot")
parser$add_argument("--followUpYears",
                    default = 5,
                    help = "Number of years of follow-up to plot")
parser$add_argument("--logFile", required=TRUE, 
                    help = "Name of (.txt) file to store logs on number of individuals plotted")
parser$add_argument("--outPlotDir", required=TRUE, 
                    help = "Path to directory to save output plots")
args <- parser$parse_args()
print(args)

# Read data ----

# UKB eids and age at diagnosis (continuous)
diag_dat <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_time_to_event_matrix.txt",
                       sep = "\t", header = T, stringsAsFactors = F)
colnames(diag_dat) <- gsub("^X", "", colnames(diag_dat))

# Disease dictionary
dictionary <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt",
                         sep = "\t", header = T, stringsAsFactors = F, 
                         quote = "", comment.char = "$")
# Subset to diagnoses of interest
DIAGS <- strsplit(args$diagnoses, ", |,|\n")[[1]]
if (DIAGS == "all") {
  DIAGS <- dictionary$phenotype
}
dictionary <- dictionary %>% filter(phenotype %in% DIAGS)

# Only keep columns of diagnosis data that are in diags
keep_ucodes <- as.character(dictionary$unique_code)
diag_dat <- diag_dat[, c("eid", keep_ucodes)]
colnames(diag_dat)[-1] <- DIAGS

BIOM <- args$biomarker
biomarker_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/indiv_qcd_data.rds")[[BIOM]] %>%
  ungroup() %>%
  mutate(eid = as.character(eid))
covar_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[[BIOM]] %>%
  ungroup() %>%
  mutate(eid = as.character(eid))
general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220131_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, comment.char = "$",
                             stringsAsFactors = F) %>%
  mutate(eid = as.character(eid))


ALL_AVAILABLE_IDS <- unique(biomarker_dat$eid)
XLIMIT <- as.numeric(args$followUpYears)

# Wrangle data and set colour palette ----

# Residualise biomarker data with respect to age (adjust for sex,
# age- and age-sq)

biomarker_dat <- biomarker_dat %>% 
  mutate(age_event_sq = age_event^2,
         sex = general_covars$sex[match(eid, general_covars$eid)])

# Return NA for missing covariates
adj_model <- lm(value ~ age_event + age_event_sq + sex, 
                    data = biomarker_dat,
                na.action = na.exclude)

modelled_dat <- biomarker_dat
modelled_dat <- modelled_dat %>% 
  mutate(adj_value = resid(adj_model),
         baseline_age = covar_dat$baseline_age[match(eid, covar_dat$eid)],
         time_from_first_mmt = age_event - baseline_age) %>%
  filter(!is.na(adj_value))

custom_col_pal <- c("#7DA15B", "#000000", "#E41A1C")
names(custom_col_pal) <- c("never_diag", "pre_diag", "post_diag")

plot_dat <- lapply(DIAGS, function (d) {
  # Prep diagnosis data
  diag_info <- diag_dat[, c("eid", d)] %>%
    filter(!is.na(!!as.symbol(d))) %>%
    mutate(eid = as.character(eid))
  colnames(diag_info) <- c("eid", "age_at_diag")
  
  # Merge in disease status (0 before diagnosis and 1 after)
  res <- modelled_dat %>% 
    mutate(age_at_diag = diag_info$age_at_diag[match(eid, diag_info$eid)],
           disease_status = factor(ifelse(is.na(age_at_diag), "never_diag",
                                          ifelse(age_event < age_at_diag, 
                                          "pre_diag", "post_diag")), 
                                   levels = c("never_diag", "pre_diag", "post_diag")))
  
  # Log to file
  sink(args$logFile, append = T)
  cat(paste0("\t", "** DIAGNOSIS ** ", d, "\n",
             "\t", "Observations retained for plot: "))
  print(table(res$disease_status))
  cat("\n")
  sink()
  return (res)
})
names(plot_dat) <- DIAGS

# Functions to plot random individuals ----

## Get random ids with and without diagnosis ----

getRandIDS <- function (diag, n_sample = 25) {
  
  unique_diag_ids <- plot_dat[[diag]] %>% 
    filter(disease_status != "never_diag") %>%
    distinct(eid)
  unique_diag_ids <- unique_diag_ids$eid
  unique_nondiag_ids <- ALL_AVAILABLE_IDS[!ALL_AVAILABLE_IDS %in% 
                                            unique_diag_ids]
  sampled_diag_ids <- sample(unique_diag_ids, min(n_sample, 
                                                  length(unique_diag_ids)), 
                             replace = F)
  sampled_nondiag_ids <- sample(unique_nondiag_ids, min(n_sample, 
                                                        length(unique_nondiag_ids)), 
                                replace = F)
  
  return (list(diag = sampled_diag_ids,
               nondiag = sampled_nondiag_ids))
}

## Functions to create plots ----

# One sub-plot per individual, faceted by id
plotIndivObsDat <- function (diag, ids_to_plot) {
  # Get biomarker data to plot and subset to relevant ids
  if (length(ids_to_plot) > 0) {
    to_plot <- plot_dat[[diag]] %>% filter(eid %in% ids_to_plot) %>%
      mutate(eid_f = factor(as.character(eid)))
    
    res <- ggplot(to_plot, aes(x = time_from_first_mmt, y = adj_value,
                               group = eid, colour = disease_status)) +
      facet_wrap(.~eid_f, ncol = 5) +
      geom_point() +
      geom_line() +
      scale_color_manual(values = custom_col_pal) +
      labs(x = "Time from first measurement (years)", 
           y = paste0("Age-adjusted ", BIOM), 
           title = diag) +
      theme(title = element_text(size = 8),
            axis.title = element_text(size = 8))
  } else res <- NULL
  return (res)
}

# Sample up to 20 individuals and plot them all 
plotManyObsDat <- function (diag) {
  # Get up to 20 random IDs (10 diagnosed, 10 undiagnosed)
  ids_to_plot <- getRandIDS(diag, 10)
  # Choose up to 5 diagnosed ids to highlight
  highlight_ids <- sample(ids_to_plot$diag, 
                          min(5, length(ids_to_plot$diag)),
                          replace = F)
  
  # Get data
  to_plot <- plot_dat[[diag]] %>% filter(eid %in% unlist(ids_to_plot)) %>%
    mutate(highlight = eid %in% highlight_ids)
  
  # Plot
  res_plot <- ggplot(to_plot %>% filter(!highlight), 
                     aes(x = time_from_first_mmt, 
                         y = adj_value, group = eid,
                         colour = disease_status)) +
    geom_point(alpha = 0.1) +
    geom_line(alpha = 0.1) +
    geom_point(data = to_plot %>% filter(highlight),
               alpha = 1) +
    geom_line(data = to_plot %>% filter(highlight),
              alpha = 1) +
    scale_color_manual(values = custom_col_pal) +
    labs(x = "Time from first measurement (years)", 
         y = paste0("Age-adjusted ", BIOM), 
         title = diag) +
    theme(title = element_text(size = 8),
          axis.title = element_text(size = 8))
  return (res_plot)
}

# Functions to plot group mean trajectories pre-diagnosis and post-diagnosis ----

# Create baselined variable (baseline pre-diagnosis, baseline post-diagnosis)
# with time from baseline

getBaselinedDat <- function (diag) {
  df <- plot_dat[[diag]] %>% 
    # Arrange by age at event 
    group_by(eid) %>%
    arrange(age_event, .by_group = T) %>%
    # Add different baseline for undiagnosed, pre- and post-diag
    ungroup() %>% group_by(eid, disease_status) %>%
    mutate(baseline = first(age_event),
           time_from_baseline = age_event - baseline)
  return (df)
}

# Get binned summaries by time spent in disease status

# Given a dataframe with times, add time bins cut by specified interval length,
# summarise values across individuals in each time bin
# and get the rolling average to plot

summTimeBin <- function (df, interval_yrs = 0.25, roll_years = 2) {
  # Cut times into 3-month (0.25 yr) bins
  cut_points <- seq(0, 80, by = interval_yrs)
  cut_labels <- seq(0, 80 - interval_yrs, by = interval_yrs)
  df$time_bin <- cut(df$time_from_baseline, 
                     breaks = cut_points, 
                     include.lowest = T,
                     labels = cut_labels)
  res <- df %>% 
    group_by(disease_status, time_bin) %>% 
    summarise(mean_value = mean(adj_value),
              sd_value = sd(adj_value), 
              n = n()) %>%
    mutate(lci_mean = mean_value - 1.96*(sd_value/sqrt(n)),
           uci_mean = mean_value + 1.96*(sd_value/sqrt(n))) 
  
  fn_roll <- function (x) { rollapply(x, roll_years/interval_yrs, 
                                      mean, fill = NA) }
  res <- res %>% ungroup() %>% 
    group_by(disease_status) %>% 
    mutate(time_bin = as.numeric(as.character(time_bin))) %>%
    arrange(time_bin, .by_group = T) %>%
    mutate(across(c(mean_value, lci_mean, uci_mean),
                  fn_roll))
  return (res)
}

## Function to create plots ----

plotMeanObsDat <- function (diag) {
  # Get data
  to_plot <- summTimeBin(getBaselinedDat(diag), 
                         interval_yrs = 0.25, roll_years = 2)
  
  res_plot <- ggplot(to_plot,
                     aes(x = time_bin, y = mean_value,
                         fill = disease_status,
                         colour = disease_status)) +
    geom_line() +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.2) +
    scale_color_manual(values = custom_col_pal) +
    scale_fill_manual(values = custom_col_pal) +
    scale_x_continuous(limits = c(0, XLIMIT)) +
    labs(x = "Time from first measurement (years)", 
         y = paste0("Mean age-adjusted ", BIOM), 
         title = diag) +
    theme(title = element_text(size = 8),
          axis.title = element_text(size = 8))
  return (res_plot)
}

# Run through each biomarker and arrange plots ----

# Random sample of individual ids
indiv_plots <- lapply(DIAGS, function (d) {
  plotIndivObsDat(d, getRandIDS(d, 5)$diag)
})

# Random sample (hundreds of ids)
many_plots <- lapply(DIAGS, function (d) {
  plotManyObsDat(d)
})

# Mean observed data pre- and post-diagnosis
mean_plots <- lapply(DIAGS, function (d) {
  plotMeanObsDat(d)
})


pdf(paste0(args$outPlotDir, BIOM, "_trajectories.pdf"), onefile = T)
ggarrange(plotlist = indiv_plots, 
          ncol = 1, nrow = 3, common.legend = T)
ggarrange(plotlist = many_plots, 
          ncol = 2, nrow = 3, common.legend = T)
ggarrange(plotlist = mean_plots, 
          ncol = 2, nrow = 3, common.legend = T)
dev.off()
