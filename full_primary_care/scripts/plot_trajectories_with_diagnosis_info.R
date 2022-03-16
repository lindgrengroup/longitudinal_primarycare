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
parser$add_argument("--idFile", required=TRUE,
                    help = "Two-column whitespace separated .txt file with id and age at diagnosis; first line is a header")
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

# UKB eids and age at diagnosis (continuous)
diag_dat <- read.table(args$idFile, header = T, stringsAsFactors = F)
colnames(diag_dat) <- c("eid", "age_at_diag")
# Clean to ensure that there are numeric values in age at diagnosis
diag_dat <- diag_dat %>% 
  filter(!is.na(age_at_diag)) %>%
  mutate(age_at_diag = as.numeric(age_at_diag))

sink(args$logFile, append = T)
cat(paste0("Number of ids provided: ", nrow(diag_dat), "\n"))
sink()

# Biomarkers longitudinal GP data
BIOMS <- strsplit(args$biomarkers, ", |,|\n|\\s")[[1]]
if (BIOMS == "all") {
  BIOMS <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/code_lists/qcd_traits_available.txt",
                      header = F, stringsAsFactors = F)$V1
}
biomarker_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/indiv_qcd_data.rds")[BIOMS]

# Wrangle data and set colour palette ----

custom_col_pal <- c("#000000", "#E41A1C")
names(custom_col_pal) <- c("pre_diag", "post_diag")

# Subset GP data to ids that are in the provided id file and 
# merge with age at diagnosis
plot_dat <- lapply(BIOMS, function (bm) {
  df <- biomarker_dat[[bm]] %>% ungroup()
  res <- df %>% filter(eid %in% diag_dat$eid) %>%
    # Merge in age at diagnosis
    # and disease status (0 before diagnosis and 1 after)
    mutate(age_at_diag = diag_dat$age_at_diag[match(eid, diag_dat$eid)],
           disease_status = factor(ifelse(age_event < age_at_diag, 
                                          "pre_diag", "post_diag"), 
                                   levels = c("pre_diag", "post_diag")))
  
  # Log to file
  sink(args$logFile, append = T)
  cat(paste0("\t", "** BIOMARKER **", bm, "\n",
             "\t", "Number of ids retained for plot: ", 
             length(unique(res$eid)), "\n"))
  sink()
  
  return (res)
})
names(plot_dat) <- BIOMS

# Functions to plot random individuals ----

## Get random ids within each class ----

getRandIDS <- function (bm, n_sample = 25) {
  unique_ids <- unique(plot_dat[[bm]]$eid)
  sampled_ids <- sample(unique_ids, min(n_sample, length(unique_ids)), 
                        replace = F)
  return (sampled_ids)
}

## Functions to create plots ----

# One sub-plot per individual, faceted by id
plotIndivObsDat <- function (bm, ids_to_plot) {
  # Get biomarker data to plot and subset to relevant ids
  to_plot <- plot_dat[[bm]] %>% filter(eid %in% ids_to_plot) %>%
    mutate(eid_f = factor(as.character(eid)))
  
  res <- ggplot(to_plot, aes(x = age_event, y = value,
                             group = eid, colour = disease_status)) +
    facet_wrap(.~eid_f, ncol = 5, scales = "free_x") +
    geom_point() +
    geom_line() +
    scale_color_manual(values = custom_col_pal) +
    labs(x = "Age (years)", y = bm)
  return (res)
}

# Sample up to 20 individuals and plot them all 
plotManyObsDat <- function (bm) {
  # Get up to 100 random IDs 
  ids_to_plot <- getRandIDS(bm, 20)
  # Choose up to 5 random ids to highlight
  highlight_ids <- sample(ids_to_plot, min(5, length(ids_to_plot)),
                          replace = F)
  
  # Get data
  to_plot <- plot_dat[[bm]] %>% filter(eid %in% ids_to_plot) %>%
    mutate(highlight = eid %in% highlight_ids)
  
  # Plot
  res_plot <- ggplot(to_plot %>% filter(!highlight), 
                     aes(x = age_event, y = value, group = eid,
                         colour = disease_status)) +
    geom_point(alpha = 0.1) +
    geom_line(alpha = 0.1) +
    geom_point(data = to_plot %>% filter(highlight),
               alpha = 1) +
    geom_line(data = to_plot %>% filter(highlight),
              alpha = 1) +
    scale_color_manual(values = custom_col_pal) +
    labs(x = "Age (years)",
         y = bm, title = bm)
  return (res_plot)
}

# Functions to plot group mean trajectories pre-diagnosis and post-diagnosis ----

# Create baselined variable (baseline pre-diagnosis, baseline post-diagnosis)
# with time from baseline

getBaselinedDat <- function (bm) {
  df <- plot_dat[[bm]] %>% 
    # Arrange by age at event 
    group_by(eid) %>%
    arrange(age_event, .by_group = T) %>%
    # Add different baseline for pre- and post-diag
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
    summarise(mean_value = mean(value),
              sd_value = sd(value), 
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

plotMeanObsDat <- function (bm) {
  # Get data
  to_plot <- summTimeBin(getBaselinedDat(bm), 
                         interval_yrs = 0.25, roll_years = 2)
  
  res_plot <- ggplot(to_plot,
                     aes(x = time_bin, y = mean_value,
                         fill = disease_status,
                         colour = disease_status)) +
    geom_line() +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.2) +
    scale_color_manual(values = custom_col_pal) +
    scale_fill_manual(values = custom_col_pal) +
    labs(x = "Time from first measurement (years)", 
         y = paste0("Mean and 95% CI"), 
         title = bm)
  return (res_plot)
}

# Run through each biomarker and arrange plots ----

# Random sample of individual ids
indiv_plots <- lapply(BIOMS, function (bm) {
  plotIndivObsDat(bm, getRandIDS(bm, 5))
})

# Random sample (hundreds of ids)
many_plots <- lapply(BIOMS, function (bm) {
  plotManyObsDat(bm)
})

# Mean observed data pre- and post-diagnosis
mean_plots <- lapply(BIOMS, function (bm) {
  plotMeanObsDat(bm)
})


pdf(paste0(args$outPlotDir, "trajectories.pdf"), onefile = T)
ggarrange(plotlist = indiv_plots, 
          ncol = 1, nrow = 3, common.legend = T)
ggarrange(plotlist = many_plots, 
          ncol = 2, nrow = 3, common.legend = T)
ggarrange(plotlist = mean_plots, 
          ncol = 2, nrow = 3, common.legend = T)
dev.off()
