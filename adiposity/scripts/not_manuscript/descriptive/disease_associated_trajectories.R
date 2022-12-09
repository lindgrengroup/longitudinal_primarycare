# Author: Samvida S. Venkatesh
# Date: 29/06/21

library(lme4)
library(splines)
library(tidyverse)
theme_set(theme_bw())
library(RColorBrewer)

# Read files ----

# EID x disease matrix
eid_pheno_matrix <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/eid_phenotype_matrix.txt", 
                               sep = "\t", header = T, stringsAsFactors = F)

# Disease dictionary
dictionary <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt",
                         sep = "\t", header = T, stringsAsFactors = F)
DISEASES <- dictionary$phenotype
colnames(eid_pheno_matrix)[-1] <- DISEASES

# BMI data 
raw_dat <- read.table("fitted_BMI_sex_comb.txt", sep = "\t",
                      header = T)

# BMI model
model_dat <- readRDS("best_model_BMI_sex_comb.rds")

# Model coefficients
model_coefs <- read.table("model_ranefs_BMI_sex_comb.txt", sep = "\t",
                          header = T)
model_coefs <- as.data.frame(t(t(as.matrix(model_coefs[, 1:4])) + 
                                 fixef(model_dat)))
model_coefs$eid <- rownames(model_coefs)

# Subset matrix to diseases of interest in population of interest ----

PLOT_DISEASES <- c("Anorexia and bulimia nervosa",
                   "Obesity",
                   "Diabetes",
                   "Polycystic ovarian syndrome",
                   "Primary Malignancy_Stomach",
                   "Hearing loss")

sub_mat <- eid_pheno_matrix[eid_pheno_matrix$eid %in% raw_dat$eid,
                                     c("eid", PLOT_DISEASES)]

# Plot population level trajectories in these diseases vs controls ----

## Raw data summary ----

# Bin into 1-yr age windows to plot
age_bin_cuts <- seq(20, 81, by = 1)
raw_dat$age_bin <- cut(raw_dat$age_event,
                       age_bin_cuts, include.lowest = T,
                       labels = seq(20, 80, by = 1))

# Function for binning raw data to plot
summarise_raw <- function (dis) {
  # Get relevant data ----
  tmp_df <- raw_dat
  tmp_df$disease_status <- as.logical(sub_mat[match(tmp_df$eid,
                                                    sub_mat$eid), 
                                              dis])
  if (dis == "Polycystic ovarian syndrome") {
    tmp_df <- subset(tmp_df, tmp_df$sex == "F")
  }
  
  # Summarise ----
  summ_raw <- tmp_df %>% group_by(disease_status, age_bin) %>% 
    summarise(mean = mean(value),
              count = n(),
              lwr = mean(value) - (1.96*sd(value)/sqrt(count)),
              upr = mean(value) + (1.96*sd(value)/sqrt(count)))
  summ_raw$data_type <- "observed"
  
  # Format for plotting ----
  summ_raw <- summ_raw[, c("disease_status", "age_bin",
                           "mean", "lwr", "upr", "data_type")]
  colnames(summ_raw) <- c("disease_status", "age",
                          "mean", "lwr", "upr", "data_type")
  summ_raw$age <- as.numeric(as.character(summ_raw$age))
  
  return (summ_raw)
}

## Predict data for new age basis ----

# Create basis for prediction
old_basis <- ns(raw_dat$scaled_age, 3)
AGE_PRED <- seq(20, 80, by = 1)
new_basis <- as.data.frame(predict(old_basis, scale(AGE_PRED)))
new_basis <- bind_cols(1, new_basis)

# Get model parameters for each individual
MODEL_PARS <-
  colnames(model_coefs)[-which(colnames(model_coefs) == "eid")]
for_traj <- merge(model_coefs, sub_mat, by = "eid")

# Function for prediction df 
summarise_pred <- function (dis) {
  
  # Get case / control centres ----
  coef_dat <- for_traj[, c("eid", MODEL_PARS, dis)]
  colnames(coef_dat)[ncol(coef_dat)] <- "disease_status"
  cluster_par <- coef_dat %>% group_by(disease_status) %>%
    summarise(across(all_of(MODEL_PARS),
                     list(mean = mean, sd = sd, n = length),
                     .names = "{.col}_{.fn}")) %>%
    # Pivot longer to calculate 95% CI (mean +/- 1.96*S.E.)
    pivot_longer(cols = -all_of("disease_status"),
                 names_to = c("parameter", ".value"),
                 names_pattern = "(.*)(_.*)") %>%
    mutate(lwr = `_mean` - (1.96*`_sd`/sqrt(`_n`)),
           upr = `_mean` + (1.96*`_sd`/sqrt(`_n`)))
  colnames(cluster_par) <- c("disease_status", "parameter",
                             "mean", "sd", "n", "lwr", "upr")
  
  # Matrix multiplication with basis ----
  # For matrix multiplication remove s.d. columns and pivot back wider
  # by disease status
  cluster_par <- cluster_par %>% select(-one_of(c("sd", "n"))) %>%
    pivot_wider(id_cols = parameter,
                names_from = disease_status,
                names_glue = "{disease_status}_{.value}",
                values_from = c(mean, lwr, upr))
  
  # Calculate trajectories
  trajectories <- as.matrix(new_basis) %*% as.matrix(cluster_par[, -1])
  trajectories <- as.data.frame(trajectories)
  
  # Merge in age, fixed effect of covariates ----
  trajectories$age_event <- AGE_PRED
  trajectories <- pivot_longer(trajectories,
                               cols = -all_of(c("age_event")),
                               names_to = c("disease_status", ".value"),
                               names_pattern = "(.*)(_.*)")
  colnames(trajectories) <- c("age", "disease_status",
                              "mean", "lwr", "upr")
  trajectories$disease_status <- 
    as.logical(as.numeric(trajectories$disease_status))
  # Add back the fixed effect of covariates 
  fxef <- tmp_df %>% group_by(disease_status) %>%
    summarise(mean_fxef = mean(fitted_covs))
  trajectories$fxef <- fxef$mean_fxef[match(trajectories$disease_status,
                                            fxef$disease_status)]
  
  summ_pred <- trajectories %>% mutate(full_mean = mean + fxef,
                                       full_lwr = lwr + fxef,
                                       full_upr = upr + fxef)
  summ_pred$data_type <- "predicted"
  
  # Format for plotting ----
  summ_pred <- summ_pred[, c("disease_status", "age",
                           "full_mean", "full_lwr", "full_upr", 
                           "data_type")]
  colnames(summ_pred) <- c("disease_status", "age",
                          "mean", "lwr", "upr", "data_type")
  return (summ_pred)
}

## Plot ----

pop_traj <- lapply(PLOT_DISEASES, function (dis) {
  
  # Get data
  summ_raw <- summarise_raw(dis)
  summ_pred <- summarise_pred(dis)
  plot_dat <- bind_rows(summ_raw, summ_pred)

  pop_plot <- ggplot(plot_dat, aes(x = age, y = mean,
                               colour = disease_status, 
                               fill = disease_status,
                               group = disease_status)) +
    facet_wrap(~data_type, nrow = 2) +
    geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    scale_colour_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Age (years)", y = "Mean (95% CI of mean) BMI", 
         title = dis)
  
  return (pop_plot)
})

pdf("population_trajectories.pdf", onefile = T)
print(pop_traj)
dev.off()


# Plot individual level trajectories in these diseases vs controls ----

## Sample individuals with and without disease ----

sampleIndivs <- function (dis, neach = 6) {
  tmp_df <- raw_dat 
  tmp_df$disease_status <- as.logical(sub_mat[match(tmp_df$eid, 
                                                    sub_mat$eid), 
                                              dis])
  if (dis == "Polycystic ovarian syndrome") {
    tmp_df <- subset(tmp_df, tmp_df$sex == "F")
  }
  
  # Sample ids to plot 
  sample_ids <- tmp_df %>% distinct(eid, disease_status) %>%
    group_by(disease_status) %>% sample_n(neach) %>%
    # Rename ids to reflect disease status
    mutate(plot_id = paste0(ifelse(disease_status == 1, "case",
                                   "control"), "_", row_number()))
  return (sample_ids)
}

## Raw data ----

get_raw <- function (ids_dat) {
  
  sample_raw <- subset(raw_dat, raw_dat$eid %in% ids_dat$eid)
  sample_raw <- merge(sample_raw, ids_dat, by = "eid")
  
  # Format for plotting ----
  sample_raw <- sample_raw[, c("eid", "plot_id", "disease_status",
                               "age_event", "value")]
  colnames(sample_raw) <- c("eid", "plot_id", "disease_status",
                          "age", "value")
  
  return (sample_raw)
}

## Predict data for new age basis ----

# Function for prediction df 
sample_pred <- function (ids_dat) {
  
  # Get parameters for the required ids ----
  indiv_par <- for_traj[, c("eid", MODEL_PARS)]
  indiv_par <- subset(indiv_par, indiv_par$eid %in% ids_dat$eid)
  
  # Matrix multiplication with basis ----
  # Calculate trajectories
  trajectories <- as.matrix(new_basis) %*% t(as.matrix(indiv_par[, -1]))
  trajectories <- as.data.frame(trajectories)
  colnames(trajectories) <- indiv_par$eid
  
  # Merge in age, fixed effect of covariates ----
  trajectories$age <- AGE_PRED
  trajectories <- pivot_longer(trajectories, 
                               cols = -all_of(c("age")),
                               names_to = "eid", values_to = "value")
  trajectories <- merge(trajectories, ids_dat, by = "eid")

  # Add back the fixed effect of covariates 
  trajectories$fxef <- 
    raw_dat$fitted_covs[match(trajectories$eid, raw_dat$eid)]
  
  sample_pred <- trajectories %>% mutate(full_value = value + fxef)

  # Format for plotting ----
  sample_pred <- sample_pred[, c("eid", "plot_id", "disease_status",
                               "age", "full_value")]
  colnames(sample_pred) <- c("eid", "plot_id", "disease_status",
                            "age", "value")
  return (sample_pred)
}

## Plot ----

indiv_traj <- lapply(PLOT_DISEASES, function (dis) {
  
  # Get data
  sample_ids <- sampleIndivs(dis)
  sample_raw <- get_raw(sample_ids)
  sample_pred <- sample_pred(sample_ids)
  
  indiv_plot <- ggplot(sample_pred, aes(x = age, y = value,
                                        colour = disease_status)) +
    facet_wrap(~plot_id, ncol = 3) +
    geom_line() +
    geom_point(data = sample_raw, colour = "black") +
    scale_x_continuous(limits = c(20, 80)) +
    scale_colour_brewer(palette = "Set1") +
    labs(x = "Age (years)", y = "BMI (Predicted)", 
         title = dis) 
  
  return (indiv_plot)
})

pdf("sample_individual_dat.pdf", onefile = T)
print(indiv_traj)
dev.off()


# Plot distribution of model coefficients in cases vs controls ----

## MODIFY THIS, NO NEED TO DOWNSAMPLE FOR DENSITY!!

# Get model parameters for each individual
MODEL_PARS <- 
  colnames(model_coefs)[-which(colnames(model_coefs) == "eid")]
for_traj <- merge(model_coefs, sub_mat, by = "eid")

parameter_dist <- lapply(PLOT_DISEASES, function (dis) {
  # Get case / control status
  tmp_dat <- for_traj[, c("eid", MODEL_PARS, dis)]
  colnames(tmp_dat)[ncol(tmp_dat)] <- "disease_status"
  tmp_dat$disease_status <- as.logical(tmp_dat$disease_status)
  
  tmp_dat$sex <- raw_dat$sex[match(tmp_dat$eid, raw_dat$eid)]
  if (dis == "Polycystic ovarian syndrome") {
    tmp_dat <- subset(tmp_dat, tmp_dat$sex == "F")
  }
  
  # # Downsample controls to 20,000 for plot 
  # NCTRL <- min(20000,
  #              nrow(tmp_dat) - sum(tmp_dat$disease_status))
  # choose_controls <- tmp_dat$eid[!tmp_dat$disease_status]
  # choose_controls <- sample(choose_controls, NCTRL, replace = F)
  # # Keep all cases and subset of controls
  # keep_ids <- c(tmp_dat$eid[tmp_dat$disease_status],
  #               choose_controls)
  # tmp_dat <- subset(tmp_dat, tmp_dat$eid %in% keep_ids)
  
  # Pivot longer by model parameter for plot
  tmp_dat <- pivot_longer(tmp_dat, 
                          cols = all_of(MODEL_PARS),
                          names_to = "parameter",
                          values_to = "parameter_value")
  # Median value in each parameter and disease status for intercept
  mu_dat <- tmp_dat %>% group_by(parameter, disease_status) %>%
    summarise(med_val = median(parameter_value))
  # Density plot of model coefficients
  dplot <- ggplot(tmp_dat, aes(x = parameter_value,
                               colour = disease_status,
                               fill = disease_status)) +
    facet_wrap(~parameter, ncol = 1, scales = "free") +
    geom_density(alpha = 0.2) +
    geom_vline(data = mu_dat, aes(xintercept = med_val,
                                  colour = disease_status),
               linetype = "dashed") +
    scale_fill_brewer(palette = "Set1") + 
    labs(title = dis) 
  
  return (dplot)
})
names(parameter_dist) <- PLOT_DISEASES

pdf("disease_model_parameters.pdf", onefile = T)
print(parameter_dist)
dev.off()

