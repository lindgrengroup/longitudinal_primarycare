# Author: Samvida S. Venkatesh
# Date: 02/03/21
# ADAPTED FROM ADIPOSITY CHANGE CONSORTIUM SOP

library(tidyverse)

# Read stratified raw slope and covariate data ----

raw_slopes <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/raw_slopes_and_covars.rds")

PHENOTYPES <- names(raw_slopes)
STRATA <- names(raw_slopes[[1]])
NPCs <- 21
PCs <- paste0("PC", 1:NPCs)

# Add sex-combined strata to covars files ----

raw_slopes <- lapply(PHENOTYPES, function (p) {
  all_dfs <- bind_rows(raw_slopes[[p]])
  # Split into sex- and ancestry-specific strata (as already existed)
  ss <- all_dfs %>% group_by(sex, ancestry) %>% group_split(.keep = T)
  names(ss) <- lapply(ss, function (x) paste(unique(x$ancestry), unique(x$sex),
                                             sep = "_"))
  # Also create a sex-combined stratum
  sc <- all_dfs %>% group_by(ancestry) %>% group_split(.keep = T)
  names(sc) <- lapply(sc, function (x) paste(unique(x$ancestry), "sexcomb",
                                             sep = "_"))
  res <- c(ss, sc)
  return (res)
})
names(raw_slopes) <- PHENOTYPES
STRATA <- names(raw_slopes[[1]])

# Build final models ----

formula_strings <- list()
# Model covariates for BMI and WHR: 
# baseline age and age-sq, baseline trait, 
# # FU years, # FU measures,
# genetic PCs, (sex)
formula_strings$mLIN <- paste0("linear_slope ~ baseline_age + age_sq + baseline_trait + FUyrs + FU_n + ", 
                             paste(PCs, collapse = " + "))
formula_strings$mLINsex <- paste0("linear_slope ~ baseline_age + age_sq + sex + baseline_trait + FUyrs + FU_n + ", 
                                paste(PCs, collapse = " + "))

formula_strings$mSPL <- paste0("spline_slope ~ baseline_age + age_sq + baseline_trait + FUyrs + FU_n + ", 
                               paste(PCs, collapse = " + "))
formula_strings$mSPLsex <- paste0("spline_slope ~ baseline_age + age_sq + sex + baseline_trait + FUyrs + FU_n + ", 
                                  paste(PCs, collapse = " + "))

# Run adjustment models and get residuals ----

adj_slopes_and_covars <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    df <- raw_slopes[[p]][[s]]
    # Remove outliers for linear and spline slopes respectively
    linear_df <- subset(df, !df$linear_slope_outlier_flag)
    spline_df <- subset(df, !df$spline_slope_outlier_flag)
    
    # Run adjustment models
      if (grepl("sexcomb", s)) { # Sex-combined model
        # Run linear and spline slope adjustments
        lin_model <- lm(formula(formula_strings$mLINsex), linear_df)
        spl_model <- lm(formula(formula_strings$mSPLsex), spline_df)
        } else { # Sex-specific model
          lin_model <- lm(formula(formula_strings$mLIN), linear_df)
          spl_model <- lm(formula(formula_strings$mSPL), spline_df)
        } 
    
    # Add residuals to covariates dataframe
    linear_df$adj_linear_slope <- lin_model$residuals 
    spline_df$adj_spline_slope <- spl_model$residuals
    
    # Combine linear and pspline adjustments with original covariates df
    res <- merge(df, linear_df[, c("eid", "adj_linear_slope")],
                 by = "eid", all.x = T)
    res <- merge(res, spline_df[, c("eid", "adj_spline_slope")],
                 by = "eid", all.x = T)
    
    return (res)
  })
  names(res) <- STRATA
  return (res)
})
names(adj_slopes_and_covars) <- PHENOTYPES

saveRDS(adj_slopes_and_covars, "/well/lindgren/UKBIOBANK/samvida/adiposity/adj_slopes_and_covars.rds")
