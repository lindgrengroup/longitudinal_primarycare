Scripts in this folder:

Scripts to perform trajectory GWAS as outlined in https://github.com/OpenMendel/TrajGWAS.jl (Ko et al. 2022). Returns results for beta (mean), tau (variance), and joint effects of beta and tau for each SNP on trait. 
1. **1_sample_qc.R** - Prepare data and covariates in long format for TrajGWAS.
2. **2_fit_null_model.jl** - Use TrajGWAS Julia package to fit mixed model with fixed and random effects (including covariates and random intercept and random time-slope within each individual) but without any SNPs - null model. 
3. **3_spatest.jl** - Perform GWAS on each chromosome with the saddle point approximation implemented in Julia TrajGWAS pacakge.
4. **submit_4_plot_results.sh** and **4_plot_gwas_results.R** - Filter SNPs and plot QQ-plots and Manhattan plots.
