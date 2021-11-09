# Author: Samvida S. Venkatesh
# Date: 06/10/21

library(lme4)
library(splines)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(grid)
theme_set(theme_bw())

set.seed(061021)

# Parse in phenotype argument
args <- commandArgs(trailingOnly = T)
PHENO <- args[1]

SEX_STRATA <- c("F", "M", "sex_comb")

log_file <- paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/log_files/regression_splines_",
                   PHENO, ".txt")

# Read files ----

dat <- 
  readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/TRAINING_SET_adiposity_split_sex.rds")[[PHENO]]
SEX_STRATA <- names(dat)

PCs <- paste0("PC", 1:21)
COVARS <- c("baseline_age", "age_sq")
# There is only one data provider for WHR so this can't be in the WHR models
if (PHENO != "WHR") {
  COVARS <- c(COVARS, "data_provider")
} 

# Wrangle data ----

# Get time from baseline measurement, time^2, and time^3
model_dat <- lapply(dat, function (df) {
  res <- df %>% mutate(t = age_event - baseline_age,
                       t2 = t^2, 
                       t3 = t^3)
  return (res)
})

# Model formula ----

# Function for natural cubic spline effect of age, with 3 degrees of freedom
# in the fixed effect
# and allow for random effects in intercept, time, time^2, and time^3
# but pare down the effects if there are too many

get_formula <- function (n_re, sx) {
  # Adjust for baseline covariates
  if (sx == "sex_comb") { mod_covars <- c(COVARS, "sex") } 
  else { mod_covars <- COVARS }
  
  covars_form <- paste0("value ~ ", 
                        paste0(mod_covars, collapse = " + "), " + ",
                        paste0(PCs, collapse = " + "))
  
  # Add spline effects to formula
  fe_spline <- "splines::ns(t, 3)"
  if (n_re == 3) re_term <- "(1 + t + t2 + t3 | eid)"
  else if (n_re == 2) re_term <- "(1 + t + t2 | eid)"
  else if (n_re == 1) re_term <- "(1 + t | eid)"
  else if (n_re == 0) re_term <- "(1 | eid)"
  
  # Full formula
  full_form <- paste0(covars_form, " + ", 
                      fe_spline, " + ", 
                      re_term)
  return (full_form)
}

# Apply formula to each strata ----

spline_models <- lapply(SEX_STRATA, function (sx) {
  # Data 
  sub_dat <- model_dat[[sx]]
  n_re <- 3
  mod_formula <- formula(get_formula(n_re, sx))
  res <- try(lmer(mod_formula, sub_dat))
  
  while (class(res) == "try-error" & n_re > -1) {
    n_re <- n_re - 1
    mod_formula <- formula(get_formula(n_re, sx))
    res <- try(lmer(mod_formula, sub_dat))
  }
  
  sink(log_file, append = T)
  cat(paste0("Model summary in strata: ", sx, "\n"))
  print(summary(res))
  sink()
  
  return (res)
})
names(spline_models) <- SEX_STRATA

# Save models
saveRDS(spline_models, 
        paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/ncubic_splines_",
               PHENO, ".rds"))

# Save coefficient values, i.e. random effects, for individuals ----

refs <- lapply(SEX_STRATA, function (sx) {
  res <- ranef(spline_models[[sx]])$eid
  res$eid <- factor(rownames(res))
  return (res)
})
names(refs) <- SEX_STRATA
saveRDS(refs, 
        paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/random_effect_terms_",
               PHENO, ".rds"))

# Plot model predictions for sample of individuals ----

# Get sample of 25 random ids for all plots
rand_eids <- lapply(SEX_STRATA, function (sx) {
  sample_eids <- sample(unique(model_dat[[sx]]$eid),
                        25, replace = F)
  return (sample_eids)
})
names(rand_eids) <- SEX_STRATA

create_prediction_df <- function (sx, ids) {
  sub_dat <- subset(model_dat[[sx]], model_dat[[sx]]$eid %in% ids)
  # Calculate maximum time-point to predict tofor each eid
  sub_dat <- sub_dat %>% group_by(eid) %>% 
    mutate(max_t = max(t))
  # Get covariates
  get_cols <- c("eid", COVARS, PCs, "max_t")
  if (sx == "sex_comb") {
    get_cols <- c(get_cols, "sex")
  }
  
  # Create new data to predict from
  new_data <- sub_dat[, get_cols] %>% distinct() 
  # Timepoints to extend to 
  t <- unlist(lapply(new_data$max_t, function (x) { 
    seq(0, x, length.out = 100) }))
  tmp <- data.frame(eid = rep(new_data$eid, each = 100),
                       t = t)
  new_data <- merge(tmp, new_data, by = "eid") %>% 
    mutate(t2 = t^2, t3 = t^3)
  
  fitted_results <- as.data.frame(predict(spline_models[[sx]], 
                                          newdata = new_data))
  colnames(fitted_results) <- "fit"
  pred_df <- bind_cols(new_data, fitted_results)
  return (pred_df)
}

# Random sample of individuals
random_pred_plots <- lapply(SEX_STRATA, function (sx) {
  sub_raw <- subset(model_dat[[sx]], model_dat[[sx]]$eid %in% rand_eids[[sx]])
  
  pred_df <- create_prediction_df(sx, rand_eids[[sx]])
  
  res <- ggplot(pred_df, aes(x = t, y = fit, group = eid)) +
    facet_wrap(~eid, ncol = 5) +
    geom_line() +
    geom_point(data = sub_raw, aes(y = value)) +
    labs(x = "Time from baseline measurement (years)",
         y = PHENO, 
         title = paste0("Randomly selected individuals in sex strata: ", sx))
  return (res)
})

# Interpret random effect terms by plotting trajectories for individuals
# with highest and lowest values for each random effect term

plot_extreme_ids <- function (sx, term, top = T) {
  term_string <- ifelse(top, paste0("high ", term),
                        paste0("low ", term))
  
  tmp <- order(refs[[sx]][, term], decreasing = top)[1:5]
  extr_eids <- refs[[sx]]$eid[tmp]
  extr_raw <- subset(model_dat[[sx]], model_dat[[sx]]$eid %in% extr_eids)
  extr_raw$type <- term_string
  
  # ALso extract a random sample as background for comparison
  ran_raw <- subset(model_dat[[sx]], model_dat[[sx]]$eid %in% rand_eids[[sx]])
  ran_raw$type <- "random"
  
  # Get predictions for all individuals in either random set or extreme set
  pred_df <- create_prediction_df(sx, 
                                  c(as.character(extr_eids), 
                                    as.character(rand_eids[[sx]])))
  pred_df$type <- ifelse(pred_df$eid %in% extr_eids, 
                         term_string, "random")
  
  # Plot
  res <- ggplot(subset(pred_df, pred_df$type != "random"), 
                aes(x = t, y = fit, group = eid, color = eid)) +
    geom_line(data = subset(pred_df, pred_df$type == "random"),
              color = "grey") +
    geom_line() +
    geom_point(data = extr_raw, aes(y = value)) +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "none") +
    labs(x = "Time from baseline measurement (years)",
         y = PHENO, 
         title = paste0(term_string, " in sex strata: ", sx))
  return (res)
}

pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/ncubic_splines_predictions_",
           PHENO, ".pdf"))

print(random_pred_plots)

lapply(SEX_STRATA, function (sx) {
  # Go through each term (intercept, time, time^2, etc.) and plot
  # both top and bottom ids 
  terms_to_plot <- colnames(refs[[sx]])[-which(colnames(refs[[sx]]) == "eid")]
  lapply(terms_to_plot, function (tm) {
    top_plot <- plot_extreme_ids(sx, tm, top = T)
    bottom_plot <- plot_extreme_ids(sx, tm, top = F)
    g <- grid.arrange(top_plot, bottom_plot, nrow = 2)
    grid.draw(g)
    return ()
  })
  return ()
})

dev.off()


