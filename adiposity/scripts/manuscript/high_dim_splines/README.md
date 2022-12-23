Scripts in this folder:

## MAIN ANALYSES

**0_prep_data.R** - Create outcomes for spline modelling of obesity traits over time by: (1) setting maximum follow-up length (i.e. 7500 days), (2) adjusting outcome for covariates (baseline age, baseline age-squared, year of birth, sex, and data provider) in sex-specific and sex-combined strata, and (3) standardising by normal scaling.
1. **1_fit_hidim_splines.R** - Main functions adapted from George Nicholson. Define 100-dimensional B-spline basis incorporating large number of knots spaced across individual's observations; smooth with AR1 prior. Hyperparameters chosen by visualising fits in **supp_choose_AR1_parameters.R** below. Plot residual variances to set sigma-squared. Batch-submit these jobs per strata with **1_submit_hidim_splines.sh**.
2. **2_sample_clustering_scheme.R** - Split individuals in each strata into 80% training and 20% validation (held-out) sets and perform all cluster centroid calculations on training set alone. Across each of S = 10 iterations, calculate cluster centroids in 5,000 randomly sampled individuals by partitioning-around-medoids (PAM) clustering on a custom distance matrix (written by George Nicholson), with cluster centroids initialisation parameters L (minimum number of measurements) and M (K-tile difference in weight or BMI M-years post baseline). Check sensitivity to K, L, M parameters, random splits of training data, and comparison to validation data as described below. Get mean centroid values across iterations after re-ordering. Batch-submit these jobs per strata and at different values of K = 2:10, L = (2, 5, 10), and M = (random, 1, 2, 5, 10) with **2_submit_hidim_splines.sh**.
3. **3_soft_clustering_results.R** - Assign bootstrapped probability of individuals belonging to each cluster (for chosen clustering parameters). 

## SUPPLEMENTARY

### Choosing AR1 hyperparameters -
1. **supp_choose_AR1_parameters.R** - Test various combinations of AR1 hyperparameters, i.e. AR1_RHO, AR1_NOISE_SD, and AR1_INTERCEPT_SD to choose values for (1) above. Plot model fits to observed data for a random sample to check smoothness and regression to mean.
2. **supp_test_range_hyperparameter_fits.R** - For a random sample of 5,000 individuals in each strata, fit various combinations of AR1 hyperparameters and save the resulting models.
3. **supp_range_hyperparameter_clusters.R** - Generate cluster centroids and assign clusters to the 5,000 individuals above. 

### Chooseing clustering hyperparameters - 
1. **supp_cluster_random_splits_validation.R** - Check sensitivity of cluster allocations across random splits of the training data. Plot number of times individual is assigned to the same cluster (modal cluster) across all 10 splits.
2. **supp_plot_parameter_selection.R** - Plot silhouette scores across values of K, L, and M hyperparameters described above.
3. **supp_table_training_vs_validation.R** - Generate descriptive characteristics (summary statistics of number of follow-up measures, length of follow-up, sex stratification, etc.) within each cluster in training and validation sets.


### Other -
1. **supp_plot_model_fits.R** - Plot model fits to observed data.
2. **supp_plot_BMI_vs_weight.R** - Calculate correlations between BLUPS from modelling of BMI and weight in each strata.
