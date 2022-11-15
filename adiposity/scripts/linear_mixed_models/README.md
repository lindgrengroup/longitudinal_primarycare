Scripts in this folder:

## MAIN ANALYSES

1. **1_apply_lmms.R** - Build linear mixed effects models for effect of time (from baseline measurement) on trait, allowing for fixed and random intercepts and slopes; adjust for covariates (baseline age, age-squared, year-of-birth, sex, and data provider). Both sex-specific and sex-combined analyses. Calculate best linear unbiased predictor (BLUP) for intercept and slope as the fixed + random effect for each individual. 

## SUPPLEMENTARY (in manuscript)

1. **supp_plot_lmm_predictions.R** - Plot model fits to observed data; separate these by covariate values to check for any systematic biases in modelling.


## SUPPLEMENTARY (not in manuscript)

1. **supp_plot_BLUP_distributions.R** - Plots to look at distributions of BLUPs from models in (1), and distribution of BLUPs within covariate groups that were adjusted for in (1) - informs the decision to adjust for covariates in GWAS.
2. **supp_plot_BLUP_interactions.R** - Plots to look at the relationship between b0 and b1, i.e. BLUP for intercept and slope from models in (1) above - informs the decision to adjust for b0 in GWAS.
