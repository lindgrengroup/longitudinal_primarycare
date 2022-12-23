# Author: Samvida S. Venkatesh
# Date: 17/05/22

library(tidyverse)

# Load data ----

mainpath <- "" # REDACTED
gen_resources_path <- "" # REDACTED
hidim_mods_path <- "" # REDACTED

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

dat <- readRDS(paste0(mainpath, "/indiv_qcd_data.rds"))[PHENOTYPES]
covars <- readRDS(paste0(mainpath, "/covariates.rds"))[PHENOTYPES]

general_covars <- read.table(paste0(gen_resources_path, "/220504_QCd_demographic_covariates.txt"),
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

MOD_COVARS <- c("baseline_age", "age_sq") # Add sex as covariate for sex-combined analyses
ADD_COVARS <- c("year_of_birth") # Add data provider as covariate if there is more than one data provider

MAX_N_DAYS <- 7500 # Number of days post baseline to be included (~20 years)

# Wrangle data (add covariates) and baseline variables ----

all_cols_gather <- c("eid", MOD_COVARS, "sex",
                     ADD_COVARS, "data_provider", "baseline_trait")

full_dat <- lapply(PHENOTYPES, function (p) {
  res <- dat[[p]] %>%
    mutate(eid = as.character(eid)) %>% 
    group_by(eid) %>%
    arrange(event_dt, .by_group = T) %>%
    mutate(t_diff = as.numeric(event_dt - first(event_dt)) + 1,
           age_t1 = first(age_event),
           value_t1 = first(value)) %>% 
    # Subset to only selected interval of days
    filter(t_diff <= MAX_N_DAYS)
  
  res <- left_join(res,
                   covars[[p]] %>% 
                     mutate(eid = as.character(eid)) %>% 
                     select(any_of(all_cols_gather)),
                   by = "eid")
  res <- left_join(res, general_covars %>% select(any_of(all_cols_gather)),
                   by = "eid")
  
  res <- res[complete.cases(res %>% select(any_of(all_cols_gather))), ]
  
  return (res)
})
names(full_dat) <- PHENOTYPES

# Residualise observed values by adjusting for confounders, and standardise (normal) ----

model_dat <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    
    # Get data to model
    df <- full_dat[[p]]
    if (sx != "sex_comb") df <- df %>% filter(sex == sx)
    
    # Get covariate sets
    minimal_adj_covars <- MOD_COVARS
    full_adj_covars <- c(MOD_COVARS, ADD_COVARS)
    if (sx == "sex_comb") {
      minimal_adj_covars <- c(minimal_adj_covars, "sex")
      full_adj_covars <- c(full_adj_covars, "sex")
    }
    if (length(unique(df$data_provider)) > 1) 
      full_adj_covars <- c(full_adj_covars, "data_provider")
    
    # Minimally adjusted model
    minmod <- lm(formula(paste0("value ~ ", 
                             paste0(minimal_adj_covars, collapse = " + "))),
                 data = df)
    # Mean and variance of residuals in minimally adjusted model
    resid_min <- residuals(minmod)
    mu_min <- mean(resid_min)
    var_min <- var(resid_min)
    
    # Full model
    fullmod <- lm(formula(paste0("value ~ ", 
                                paste0(full_adj_covars, collapse = " + "))),
                 data = df)
    # Mean and variance of residuals in fully adjusted model
    resid_full <- residuals(fullmod)
    mu_full <- mean(resid_full)
    var_full <- var(resid_full)
    
    return (list(minmod = minmod, mu_min = mu_min, var_min = var_min,
                 fullmod = fullmod, mu_full = mu_full, var_full = var_full))
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(model_dat) <- PHENOTYPES

saveRDS(model_dat, 
        paste0(hidim_mods_path, "/data/models_for_refitting.rds"))

to_write <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    
    # Get data 
    df <- full_dat[[p]]
    if (sx != "sex_comb") df <- df %>% filter(sex == sx)
    
    res <- df
    mod_dat <- model_dat[[p]][[sx]]
    res$value_minadj_norm <- (residuals(mod_dat$minmod) - mod_dat$mu_min)/sqrt(mod_dat$var_min)
    res$value_fulladj_norm <- (residuals(mod_dat$fullmod) - mod_dat$mu_full)/sqrt(mod_dat$var_full)
    
    to_save <- res[, c("eid", "t_diff", 
                       "value", "value_minadj_norm", "value_fulladj_norm", 
                       "age_t1", "value_t1")]
    return (to_save)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(to_write) <- PHENOTYPES

saveRDS(to_write, 
        paste0(hidim_mods_path, "/data/dat_to_model_standardised.rds"))
