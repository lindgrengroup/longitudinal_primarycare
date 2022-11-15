Scripts in this folder:

## MAIN ANALYSES

**0_prep_data.R** - Create outcomes for spline modelling of obesity traits over time by: (1) setting maximum follow-up length (i.e. 7500 days), (2)adjusting outcome for covariates (baseline age, baseline age-squared, year of birth, sex, and data provider) in sex-specific and sex-combined strata, and (3) standardising by normal scaling.
1. **1_fit_hidim_splines.R** - Main functions adapted from George Nicholson. Define 100-dimensional B-spline basis incorporating large number of knots spaced across individual's observations; smooth with AR1 prior. Hyperparameters chosen by visualising fits in *supp_choose_AR1_parameters.R* below. Plot residual variances to set sigma-squared.
2. 

## SUPPLEMENTARY (in manuscript)

1. **supp_choose_AR1_parameters.R** - Test various combinations of AR1 hyperparameters, i.e. AR1_RHO, AR1_NOISE_SD, and AR1_INTERCEPT_SD to choose values for (1) above. Plot model fits to observed data for a random sample to check smoothness and regression to mean.
2. **supp_plot_model_fits.R** - Plot model fits to observed data.


## SUPPLEMENTARY (not in manuscript)

1. 