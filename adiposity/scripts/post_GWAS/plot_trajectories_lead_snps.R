# Author: Samvida S. Venkatesh
# Date: 28/02/2022

library(lme4)
library(splines)
library(zoo)
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

# Read files ----

args <- commandArgs(trailingOnly = T)
PARAMETER <- args[1]
STRATA <- args[2]

p <- gsub("_.*", "", STRATA)
sx <- gsub(paste0(p, "_"), "", STRATA)

# Variants of interest
lead_snps <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/post_GWAS/lead_snps/",
                               STRATA, "_", PARAMETER, "_final.lead_snps.txt"),
                        sep = "\t", header = T, stringsAsFactors = F)
VARIDS <- unique(lead_snps$rsid)

## FIRST TIME THIS SCRIPT IS RUN, COMBINE ALL GENETIC DATA INTO A SINGLE FILE
## COMMENT OUT AFTERWARDS

# gen_dat <- lapply(VARIDS, function (varid) {
#   res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/sample_variant_counts/",
#                            PARAMETER, "_", varid, "_dosages.txt"), 
#                     sep = " ", header = T, stringsAsFactors = F)
#   # Remove first row, which contains info on type of column and column 2 (ID repeat)
#   res <- res[-1, -2]
#   colnames(res) <- c("eid", "missing", "sex", varid)
#   return (res)
# })
# names(gen_dat) <- VARIDS
# 
# # Merge genetic data files
# gen_dat <- gen_dat %>% reduce(full_join, by = c("eid", "missing", "sex"))
# 
# # Write table
# write.table(gen_dat, 
#             paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/sample_variant_counts/",
#             PARAMETER, "_ALL_variants_dosages.txt"),
#             sep = "\t", row.names = F, quote = F)

# Dosage file
gen_dat <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/sample_variant_counts/lmm_slopes_adj_baseline_ALL_variants_dosages.txt",
                      sep = "\t", header = T, stringsAsFactors = F)

slope_models <- readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/lmm/lmm_",
                               p, ".rds"))[[sx]]

covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/covariates.rds")[[p]]

PCs <- paste0("PC", 1:21)
COVARS <- c("baseline_age", "age_sq", "FUyrs", "data_provider")

# Wrangle data ----

# Add in covariates to raw GP data and subset to ids that actually made it
# to the models
# model_dat <- merge(dat, covars, by = "eid")
# if (sx != "sex_comb") model_dat <- model_dat %>% filter(sex == sx)
# model_dat <- model_dat %>% mutate(t = age_event - baseline_age) 

model_dat <- getData(slope_models)

# Subset dosage data to variants of interest in samples of interest
# and pivot wider
gen_dat <- gen_dat %>% select(all_of(c("eid", VARIDS))) %>% 
  filter(eid %in% model_dat$eid) %>%
  pivot_longer(cols = all_of(VARIDS),
               names_to = "variant", values_to = "dosage")

# Convert dosages to 0/1/2 genotypes based on threshold
DOSAGE_THRESHOLD_0 <- 0.5
DOSAGE_THRESHOLD_2 <- 1.5
gen_dat <- gen_dat %>% 
  mutate(genotype = ifelse(dosage < DOSAGE_THRESHOLD_0, 0,
                           ifelse(dosage > DOSAGE_THRESHOLD_2, 2, 1)))

# Count number of individuals with each variant and genotype (0,1,2)
# to make sampling quicker later
gen_dat <- gen_dat %>% group_by(variant, genotype) %>% mutate(n_group = n())

# Plot for random individuals within each dosage group (multipage) ----
## Get random IDs within each group (variant, genotype) ----

get_rand_eids <- function (varname, n_each = 5) {
  # Return df of id, variant, dosage
  sampled_ids <- gen_dat %>% 
    filter(variant == varname) %>%
    group_by(genotype) %>%
    sample_n(ifelse(n_group < n_each, n_group, n_each)) %>%
    select(eid, genotype)
  return (sampled_ids)
}

## Function to create predicted data for set of ids ----

create_prediction_df <- function (ids) {
  sub_dat <- subset(model_dat, 
                    model_dat$eid %in% ids)
  # Calculate maximum time-point to predict to for each eid
  sub_dat <- sub_dat %>% group_by(eid) %>% 
    mutate(max_t = ceiling(max(t)))
  
  # Create new data to predict from
  new_data <- sub_dat %>% 
    select(all_of(c("eid", "max_t", "sex", PCs, COVARS))) %>% 
    distinct(eid, data_provider, .keep_all = T) 
  # Timepoints to extend to
  ts <- lapply(1:nrow(new_data), FUN = function (i) { 
    seq(0, new_data$max_t[i], by = 0.25) })
  length_ts <- unlist(lapply(ts, function (x) length(x)))
  t <- unlist(ts)
  tmp <- data.frame(eid = rep(new_data$eid, length_ts),
                    data_provider = rep(new_data$data_provider, length_ts),
                    t = t)
  new_data <- merge(tmp, new_data, by = c("eid", "data_provider"))
  
  # Predict new values
  fitted_results <- 
    as.data.frame(predict(slope_models, newdata = new_data))
  colnames(fitted_results) <- "fit"
  pred_df <- bind_cols(new_data, fitted_results)
  
  # At each time-point, average across data providers for each individual
  # Add in data on age at predicted event for plot
  pred_df <- pred_df %>% group_by(eid, t) %>%
    summarise(fit = mean(fit)) %>%
    left_join(covars, by = "eid") %>%
    mutate(age_event = t + baseline_age) %>%
    select(eid, t, age_event, fit)
  
  return (pred_df)
}

## Function to create plots ----

plot_predictions <- function (id_df, varname) {
  # Subset to one variant at a time
  gen_var <- gen_dat %>% filter(variant == varname)
  
  # Get raw data
  raw_dat <- model_dat %>% filter(eid %in% id_df$eid)
  # Plot predicted data
  plot_dat <- create_prediction_df(id_df$eid)
  plot_dat$genotype <- gen_var$genotype[match(plot_dat$eid, gen_var$eid)]
  plot_dat$genotype <- as.factor(as.character(plot_dat$genotype))
  
  # Order the ids by genotype for plot
  id_levels <- id_df$eid[order(id_df$genotype)]
  plot_dat$eid_f <- factor(as.character(plot_dat$eid), levels = id_levels)
  raw_dat$eid_f <- factor(as.character(raw_dat$eid), levels = id_levels)
  
  res <- ggplot(plot_dat, aes(x = age_event)) +
    facet_wrap(.~eid_f, nrow = 3, scales = "free_x") +
    geom_point(data = raw_dat, aes(y = value),
               colour = "black") +
    geom_line(aes(y = fit, colour = genotype)) +
    scale_color_brewer(palette = "YlOrRd", direction = -1) +
    labs(x = "Age (years)",
         y = p,
         title = paste0("Randomly selected ", sx, ", phenotype: ", p))
  return (res)
}

# Plot mean observed data per dosage for each variant ----

get_mean_obs_dat <- function (varname) {
  sub_gen <- gen_dat %>% filter(variant == varname)
  # Observed data
  obs_dat <- model_dat %>%
    mutate(genotype = sub_gen$genotype[match(eid, sub_gen$eid)],
           age_bin = plyr::round_any(age_event, 0.25, f = floor)) 
  
  summ_dat <- obs_dat %>% 
    group_by(genotype, age_bin) %>% 
    summarise(mean_value = mean(value),
              sd_value = sd(value), 
              n = n()) %>%
    mutate(lci_mean = mean_value - 1.96*(sd_value/sqrt(n)),
           uci_mean = mean_value + 1.96*(sd_value/sqrt(n))) 
  
  summ_dat <- summ_dat %>% 
    group_by(genotype) %>%
    arrange(age_bin, .by_group = T) %>%
    mutate(interval_width = seq_along(age_bin) - 
             findInterval(age_bin - 2, age_bin),
           mean_rolled_value = rollapplyr(mean_value, interval_width,
                                          mean, partial = T, na.rm = T)) 
  
  # For plotting
  summ_dat$genotype <- as.factor(as.character(summ_dat$genotype))
  summ_dat$variant <- varname
  return (summ_dat)
}

plot_mean_obs_dat <- lapply(VARIDS, function (v) get_mean_obs_dat(v))
plot_mean_obs_dat <- bind_rows(plot_mean_obs_dat)
mean_obs_plot <- ggplot(plot_mean_obs_dat,
                    aes(x = age_bin, y = mean_value,
                        fill = genotype,
                        colour = genotype)) +
  facet_wrap(~variant, ncol = 2) +
  geom_line() +
  geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), linetype = 0,
              alpha = 0.2) +
  scale_color_brewer(palette = "YlOrRd", direction = -1) +
  scale_fill_brewer(palette = "YlOrRd", direction = -1) +
  labs(x = "Age (years)",
       y = paste0("Mean and 95% C.I. of ", p),
       title = paste0("Observed trajectories in each dosage group of phenotype: ",
                      p, " strata: ", sx))

# Re-fit models within each group and plot mean trajectory per dosage ----

refit_data_by_genotype <- function (varname) {
  sub_gen <- gen_dat %>% filter(variant == varname)
  # Create modelled data split by genotype
  to_model <- model_dat 
  to_model$genotype <- 
    sub_gen$genotype[match(to_model$eid, sub_gen$eid)]
  
  fitted_dat_within_g <- lapply(c(0, 1, 2), function (ki) {
    df <- subset(to_model, to_model$genotype == ki)
    mod_covars <- c("baseline_age", "age_sq")
    if (sx == "sex_comb") mod_covars <- c(mod_covars, "sex")
    mod_formula <- formula(paste0("value ~ t + ", 
                                  paste0(mod_covars, collapse = " + "), 
                                  " + (t | eid)"))
    mod_res <- lmer(mod_formula, data = df, REML = F)
    # Fixed effect prediction 
    newdat <- create_prediction_df(df$eid)
    newdat_sumstats <- newdat %>% group_by(t) %>% 
      summarise(mean_fit = mean(fit),
                sd_fit = sd(fit), 
                n = n()) %>%
      mutate(lci_mean = mean_fit - 1.96*(sd_fit/sqrt(n)),
             uci_mean = mean_fit + 1.96*(sd_fit/sqrt(n)),
             genotype = unique(df$genotype)) 
    return (newdat_sumstats)
  })
  plot_dat <- bind_rows(fitted_dat_within_g)
  
  # Filter plots to stop at 20 years post first measurement
  plot_dat <- plot_dat %>% filter(t <= 20)
  
  plot_dat$variant <- varname
  plot_dat$genotype <- as.factor(as.character(plot_dat$genotype))
  return (plot_dat)
}

plot_mean_modelled_dat <- lapply(VARIDS, function (v) refit_data_by_genotype(v))
plot_mean_modelled_dat <- bind_rows(plot_mean_modelled_dat)
mean_modelled_plot <- ggplot(plot_mean_modelled_dat,
                        aes(x = t, y = mean_fit,
                            fill = genotype,
                            colour = genotype)) +
  facet_wrap(~variant, ncol = 2) +
  geom_line() +
  geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), linetype = 0,
              alpha = 0.2) +
  scale_color_brewer(palette = "YlOrRd", direction = -1) +
  scale_fill_brewer(palette = "YlOrRd", direction = -1) +
  labs(x = "Time from baseline measurement (years)",
       y = paste0("Mean and 95% C.I. of ", p),
       title = paste0("Modelled trajectories in each dosage group of phenotype: ",
                      p, " strata: ", sx))

# Print plots ----

pdf(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/plots/",
           p, "_", sx, "/trajectories_lead_snps_", PARAMETER, ".pdf"),
    onefile = T)
# Go through each variant to plot random sample of individuals
lapply(VARIDS, function (v) {
  plot_predictions(get_rand_eids(v, 5), v)
})
print(mean_obs_plot)
print(mean_modelled_plot)
dev.off()
