# Author: Samvida S. Venkatesh
# Date: 08/02/2022

scripts_path <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sensitivity_longevity/scripts" # REDACTED

STRATA_NAMES <- c("BMI_F", "BMI_M", "BMI_sex_comb", 
                  "Weight_F", "Weight_M", "Weight_sex_comb")
# GWAS submission script
submission_script <- paste0(scripts_path, "/5a_run_b1_associations_BOLT.sh")

# For lmm slopes
for (s in STRATA_NAMES) {
    job_options <- paste0(
      "--export=",
      paste0(
        "STRATA=\"", s, "\""
      )
    )
    job_submission <- paste("sbatch", job_options, submission_script)
    system(job_submission)
    print(job_submission)
}

# For clusters
CLUSTERS <- c("k1", "k1_k2", "k1_k2_k3")
submission_script <- paste0(scripts_path, "/5b_run_clust_associations_BOLT.sh")

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
