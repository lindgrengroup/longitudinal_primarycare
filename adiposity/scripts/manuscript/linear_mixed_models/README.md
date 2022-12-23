Scripts in this folder:

## MAIN ANALYSES

1. **1_apply_lmms.R** - Build linear mixed effects models for effect of time (from baseline measurement) on trait, allowing for fixed and random intercepts and slopes; adjust for covariates (baseline age, age-squared, year-of-birth, sex, and data provider). Both sex-specific and sex-combined analyses. Calculate best linear unbiased predictor (BLUP) for intercept and slope as the fixed + random effect for each individual. 

## SUPPLEMENTARY

1. **supp_plot_lmm_predictions.R** - Plot model fits to observed data; separate these by covariate values to check for any systematic biases in modelling.
2. **supp_plot_BMI_vs_weight.R** - Calculate correlations between BLUPS from modelling of BMI and weight in each strata.
