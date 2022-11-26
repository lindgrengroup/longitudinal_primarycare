Scripts in this folder:

Plotting and descriptive statistics for various analuses

### Raw data - 
1. **adiposity_ascertainment.R** - Generate tables and figures for various types of correlations between mean obesity trait and covariates to describe ascertainment bias effects.
1. **adiposity_report_characteristics.R** - Generate tables and figures for descriptive characteristics of raw data, similar to above but not only for mean obesity.
2. **disease_associated_trajectories.R** - For a sample of diseases, plot raw data trajectories in individuals with and without disease.

### Modelling results - 
Change these scripts for each set of models, linear, cubic splines, etc.) - 
1. **plot_BLUP_distributions.R** - Plot randomly selected subset of BLUPs from models to check distributions as well as relationships between BLUPs and model covariates.
2. **plot_all_model_predictions.R** - Plot individual-level model predictions for randomly selected individuals, those with very few or many repeat measures, at the tails of model BLUPs, etc. for both linear and cubic spline models.

### Modelling results within covariate groups - 
Scripts in */trajectories_by_covariates/* - Plot adiposity trajectories (modelled by cubic splines) in each group of a covariate (ex. birth cohort, sex, disease chapter, etc.)