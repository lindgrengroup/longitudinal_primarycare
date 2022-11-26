Scripts in this folder:

1. **1_apply_time_splines.R** - Build cubic spline mixed effects models for effect of time on trait, allowing for fixed and random intercepts and spline effects of time. Run through models with increasing number of knots in the fixed effects (3-10) and random effects (1 - # fixed effects) to choose model with lowest BIC. Adjust for covariates (baseline age, baseline age-squared, and sex). Both sex-specific and sex-combined analyses. Calculate best linear unbiased predictor (BLUP) for each term as the fixed + random effect for each individual. 
2. **2_plot_cspline_predictions.R** - Plot model predictions for various subsets of individuals (random, high or low values of baseline trait, baseline age, covariates, etc.) 
3. **3_cBLUP_PCA_heatmaps.R** - Calculate PCs based on BLUPs from cubic spline models and plot heatmaps of BLUPs as well as PCA scree plots to visualise dimensionality reduction captured by PCs.
4. **4_*.R** - Try various clustering techniques: Gaussian mixture modelling, bootstrap number of clusters, and k-means clustering - all based on intercept and BLUPs from models.
5. **5_plot_clustering_results...** - Once clusters have been determined, plot modelled trajectories, observed trajectories, and associations between clusters and various covariates.
6. **6_remap_clusters.R** - If clusters have to be made to correspond across different strata.
