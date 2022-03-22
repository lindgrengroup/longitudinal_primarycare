# Author: Samvida S. Venkatesh
# Date: 07/03/2022

# Load an R module with Bioconductor to automatically access these packages
# Otherwise specify your own package libraries to find packages

library(argparse)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(zoo)
theme_set(theme_bw())

# Parse arguments ----

parser <- ArgumentParser()
parser$add_argument("--rsid", required=TRUE,
                    help = "rsid for which dosages have already been calculated")
parser$add_argument("--biomarkers", 
                    default = "FAI, FSH, LH, Oestradiol, Progesterone, Testosterone",
                    help = "Biomarkers to plot")
parser$add_argument("--outDir", required=TRUE, 
                    help = "Path to directory to store output logs and plots")
args <- parser$parse_args()
print(args)

logFile <- paste0(args$outDir, args$rsid, "_plotting_log.txt")

# Read data ----

# UKB eids and their allele dosages for the given rsid
id_dosage <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/pott_lead_snps/",
                               args$rsid, "_dosages.txt"),
                        header = T, stringsAsFactors = F)
# Get the relevant columns (id and dosage)
id_dosage <- id_dosage[-1, c(1,5)]
colnames(id_dosage) <- c("eid", "dosage")

# Biomarkers longitudinal GP data
BIOMS <- strsplit(args$biomarkers, ", |,|\n|\\s")[[1]]
if (BIOMS == "all") {
  BIOMS <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/code_lists/qcd_traits_available.txt",
                      header = F, stringsAsFactors = F)$V1
}
biomarker_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/indiv_qcd_data.rds")[BIOMS]
biomarker_dat <- lapply(biomarker_dat, function (df) {
  res <- df %>%
    ungroup() %>%
    mutate(eid = as.character(eid))
  return (res)
})

covar_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[BIOMS]
covar_dat <- lapply(covar_dat, function (df) {
  res <- df %>%
    ungroup() %>%
    mutate(eid = as.character(eid))
  return (res)
})

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220131_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, comment.char = "$",
                             stringsAsFactors = F) %>%
  mutate(eid = as.character(eid))

# Wrangle data ----

# Convert dosages to 0/1/2 genotypes based on threshold
DOSAGE_THRESHOLD_0 <- 0.5
DOSAGE_THRESHOLD_2 <- 1.5
GTS <- c("0", "1", "2")
id_dosage <- id_dosage %>% 
  mutate(genotype = ifelse(dosage < DOSAGE_THRESHOLD_0, "0",
                           ifelse(dosage > DOSAGE_THRESHOLD_2, "2", "1")),
         genotype = factor(as.character(genotype), 
                           levels = GTS))

sink(logFile, append = T)
cat(paste0("Number of ids provided: ", nrow(id_dosage), "\n"))
cat(paste0("Classification: ", "\n"))
print(table(id_dosage$genotype))
cat("\n")
sink()

# Subset hormone data to ids that are in the provided id file and 
# merge with classification

# Adjust biomarker values for age, age^2, and sex
plot_dat <- lapply(BIOMS, function (bm) {
  df <- biomarker_dat[[bm]] %>%
    filter(eid %in% id_dosage$eid) %>%
    mutate(age_event_sq = age_event^2,
           sex = general_covars$sex[match(eid, general_covars$eid)])
  
  # Return NA for data that is missing covariates
  adj_model <- lm(value ~ age_event + age_event_sq + sex, 
                  data = df,
                  na.action = na.exclude)
  
  modelled_dat <- df
  modelled_dat <- modelled_dat %>% 
    left_join(., covar_dat[[bm]], by = "eid") %>%
    mutate(adj_value = resid(adj_model),
           time_from_baseline = age_event - baseline_age) %>%
    filter(!is.na(adj_value))
  
  # Merge in classification
  modelled_dat$genotype <- id_dosage$genotype[match(modelled_dat$eid, id_dosage$eid)]
  modelled_dat <- modelled_dat %>% mutate(genotype = factor(as.character(genotype),
                                                         levels = GTS))
  return (modelled_dat)
})
names(plot_dat) <- BIOMS

# Subset covariate data to ids that are in the provided id file and
# merge with classification
covar_plot_dat <- lapply(BIOMS, function (bm) {
  df <- plot_dat[[bm]] %>% distinct(eid, .keep_all = T) %>%
    # Subset to relevant data and pivot longer for plot
    select(eid, genotype, 
           baseline_age, baseline_trait,
           FUyrs, FU_n) %>%
    pivot_longer(cols = all_of(c("baseline_age", "baseline_trait",
                                 "FUyrs", "FU_n")),
                 names_to = "covariate", values_to = "covar_value")
  return (df)
})
names(covar_plot_dat) <- BIOMS

# Count number of individuals in each class to make sampling quicker later
summ_dat <- lapply(BIOMS, function (bm) {
  res <- covar_plot_dat[[bm]] %>% 
    distinct(eid, genotype) %>%
    group_by(genotype) %>% 
    mutate(n_group = n())
  
  # Log to file
  sink(logFile, append = T)
  cat(paste0("\t", "** BIOMARKER **", bm, "\n",
             "\t", "Number of ids retained for plot: ", nrow(res), "\n"))
  cat(paste0("\t", "Classification: "))
  print(table(res$genotype))
  cat("\n")
  sink()
  
  return (res)
})
names(summ_dat) <- BIOMS

# Create color palette based on dosages
custom_col_pal <- brewer.pal(3, "Set1")
names(custom_col_pal) <- c("0", "1", "2")

# Function to plot proportion of individuals with each genotype ----

plotProportionGenotype <- function () {
  summ_full_ukb <- id_dosage %>% 
    group_by(genotype) %>% 
    summarise(n_group = n()) %>%
    mutate(category = "full_ukbb") %>% ungroup()
  
  summ_hormones_measured <- lapply(BIOMS, function (bm) {
    res <- summ_dat[[bm]] %>%
      distinct(n_group, genotype) %>%
      mutate(category = bm) %>% ungroup()
    return (res)
  })
  summ_hormones_measured <- bind_rows(summ_hormones_measured) 
  
  to_plot <- bind_rows(summ_full_ukb, summ_hormones_measured) %>%
    mutate(category = factor(category, levels = c(BIOMS, "full_ukbb"))) %>%
    group_by(category) %>%
    mutate(total = sum(n_group),
           freq = n_group / total)
  
  res_plot <- ggplot(to_plot, aes(x = category, y = freq,
                                  fill = genotype, colour = genotype)) +
    geom_bar(position = "dodge", stat = "identity") + 
    scale_fill_manual(values = custom_col_pal) + 
    scale_colour_manual(values = custom_col_pal) +
    labs(x = "Individuals with GP measurement", 
         y = "Fraction of population") 
  return (res_plot)
}

# Functions to plot covariate distributions within each dosage ----

plotCovarDistribution <- function (bm) {
  to_plot <- covar_plot_dat[[bm]] 
  res_plot <- ggplot(to_plot, 
                     aes(x = genotype, y = covar_value)) +
    facet_wrap(~covariate, ncol = 2, scales = "free_y") +
    geom_violin(aes(fill = genotype), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_manual(values = custom_col_pal) + 
    labs(x = paste0(args$rsid, " genotype"), 
         y = "Covariate value", 
         title = bm) +
    theme(legend.position = "none")
  return (res_plot)
}

# Functions to plot random individuals within each class ----

## Get random ids within each class ----

getRandIDS <- function (bm, n_each = 5) {
  # Return dataframe of id and class
  sampled_ids <- summ_dat[[bm]] %>%
    group_by(genotype) %>%
    sample_n(ifelse(n_group < n_each, n_group, n_each)) %>%
    select(eid, genotype)
  return (sampled_ids)
}

## Functions to create plots ----

# One sub-plot per individual, faceted by id
plotIndivObsDat <- function (bm, id_df) {
  # Get biomarker data to plot and subset to relevant ids
  to_plot <- plot_dat[[bm]] %>% filter(eid %in% id_df$eid)
  
  # Order the ids by class for plot
  id_df <- to_plot %>% distinct(eid, genotype)
  id_levels <- id_df$eid[order(id_df$genotype)]
  to_plot <- to_plot %>% mutate(eid_f = factor(as.character(eid),
                                               levels = id_levels))
  
  res <- ggplot(to_plot, aes(x = age_event, y = value,
                             group = eid)) +
    facet_wrap(.~eid_f, ncol = 5) +
    geom_point(aes(y = value), colour = "black") +
    geom_line(aes(colour = genotype)) +
    scale_color_manual(values = custom_col_pal) +
    labs(x = "Age (years)",
         y = bm, title = bm)
  return (res)
}

# Functions to plot group mean trajectories within each class ----

## Time-binned summaries ----

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
    group_by(genotype, time_bin) %>% 
    summarise(mean_value = mean(adj_value, na.rm = T))  
  
  fn_roll <- function (x) { rollapply(x, roll_years/interval_yrs, 
                                      mean, fill = NA) }
  res <- res %>% ungroup() %>% 
    group_by(genotype) %>% 
    mutate(time_bin = as.numeric(as.character(time_bin))) %>%
    arrange(time_bin, .by_group = T) %>%
    mutate(mean_value = rollapply(mean_value, roll_years/interval_yrs, 
                                  mean, fill = NA),
           sd_value = rollapply(mean_value, roll_years/interval_yrs, 
                                sd, fill = NA),
           lci_mean = mean_value - 1.96*sd_value,
           uci_mean = mean_value + 1.96*sd_value)
  return (res)
}

## Function to create plots ----

plotMeanObsDatNotBaselined <- function (bm) {
  # Get data
  to_plot <- summTimeBin(plot_dat[[bm]], 
                         interval_yrs = 0.25, roll_years = 2)
  
  res_plot <- ggplot(to_plot,
                     aes(x = time_bin, y = mean_value,
                         fill = genotype,
                         colour = genotype)) +
    geom_line() +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.2) +
    scale_color_manual(values = custom_col_pal) +
    scale_fill_manual(values = custom_col_pal) +
    labs(x = "Time from first measurement",
         y = "Mean adjusted value",
         title = bm)
  return (res_plot)
}

plotMeanObsDatBaselined <- function (bm) {
  # Get data
  baselined_dat <- plot_dat[[bm]] %>%
    group_by(eid) %>%
    arrange(time_from_baseline, .by_group = T) %>%
    mutate(first_adj_value = first(adj_value),
           adj_value = adj_value - first_adj_value)
  
  to_plot <- summTimeBin(baselined_dat, 
                         interval_yrs = 0.25, roll_years = 2)
  
  res_plot <- ggplot(to_plot,
                     aes(x = time_bin, y = mean_value,
                         fill = genotype,
                         colour = genotype)) +
    geom_line() +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.2) +
    scale_color_manual(values = custom_col_pal) +
    scale_fill_manual(values = custom_col_pal) +
    labs(x = "Time from first measurement",
         y = "Mean baselined adj. value",
         title = bm)
  return (res_plot)
}

# Run through each biomarker and print plots ----

# Covariate distribution
covar_plots <- lapply(BIOMS, function (bm) {
  plotCovarDistribution(bm)
})

# Random sample of individual ids
indiv_plots <- lapply(BIOMS, function (bm) {
  plotIndivObsDat(bm, getRandIDS(bm, 5))
})

# Mean observed data within each dosage group
mean_plots <- lapply(BIOMS, function (bm) {
  plotMeanObsDatNotBaselined(bm)
})

mean_slope_plots <- lapply(BIOMS, function (bm) {
  plotMeanObsDatBaselined(bm)
})

pdf(paste0(args$outDir, args$rsid, "_trajectories.pdf"), onefile = T)
plotProportionGenotype()
ggarrange(plotlist = covar_plots, 
          ncol = 1, nrow = 2, common.legend = T)
ggarrange(plotlist = indiv_plots, 
          ncol = 1, nrow = 2, common.legend = T)
ggarrange(plotlist = mean_plots, 
          ncol = 2, nrow = 3, common.legend = T)
ggarrange(plotlist = mean_slope_plots, 
          ncol = 2, nrow = 3, common.legend = T)
dev.off()
