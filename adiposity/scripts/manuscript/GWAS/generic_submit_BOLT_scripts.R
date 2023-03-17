# Author: Samvida S. Venkatesh
# Date: 08/02/2022
# Adapted from: Duncan Palmer (https://github.com/astheeggeggs/SAIGE_gene_munging)

# R wrapper to submit BOLT jobs
# Need to loop over all of the strata and parameter types in each strata

scripts_path <- "" # REDACTED

STRATA_NAMES <- c("BMI_F", "BMI_M", "BMI_sex_comb", 
                  "Weight_F", "Weight_M", "Weight_sex_comb")
PARAMETERS <- c("b0", "b1")
# PARAMETERS <- "b1_no_FU_adjustment"
# PARAMETERS <- "b1_unadj_b0"
# CLUSTERS <- c("k1", "k1_k2", "k1_k2_k3")

# GWAS submission script
# submission_script <- paste0(scripts_path, "/BOLT_helpers/1_perform_GWAS_BOLT.sh")
# submission_script <- paste0(scripts_path, "/softprob_linear_regression_GWAS.sh")
# submission_script <- paste0(scripts_path, "/softprob_linear_regression_GWAS_no_FU_adjustment.sh")
# submission_script <- paste0(scripts_path, "/softprob_linear_regression_GWAS_no_baseline_adjustment.sh")

# Filtering submission script
# submission_script <- paste0(scripts_path, "/BOLT_helpers/2_perform_BOLT_filtering.sh")

# Fine-mapping submission script
# submission script <- paste0(scripts_path, "/BOLT_helpers/3_perform_finemapping_BOLT.sh")

submission_script <- paste0(scripts_path, "/BOLT_helpers/1_perform_GWAS_BOLT.sh")

for (s in STRATA_NAMES) {
  for (pr in PARAMETERS) {
    job_options <- paste0(
      "--export=",
      paste0(
        "STRATA=\"", s, "\",",
        "PARAMETER=\"", pr, "\""
      )
    )
    job_submission <- paste("sbatch", job_options, submission_script)
    system(job_submission)
    print(job_submission)
  }
}

for (s in STRATA_NAMES) {
  for (k in CLUSTERS) {
    job_options <- paste0(
      "--export=",
      paste0(
        "STRATA=\"", s, "\",",
        "CLUST=\"", k, "\""
      )
    )
    job_submission <- paste("sbatch", job_options, submission_script)
    system(job_submission)
    print(job_submission)
  }
}
