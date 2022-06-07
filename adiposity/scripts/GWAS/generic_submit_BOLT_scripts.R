# Author: Samvida S. Venkatesh
# Date: 08/02/2022
# Adapted from: Duncan Palmer (https://github.com/astheeggeggs/SAIGE_gene_munging)

# R wrapper to submit BOLT jobs
# Need to loop over all of the strata and parameter types in each strata

STRATA_NAMES <- c("BMI_F", "BMI_M", "BMI_sex_comb", 
                  "Weight_F", "Weight_M", "Weight_sex_comb")
PARAMETERS <- c("lmm_intercepts", "lmm_slopes_adj_int", 
                "lmm_slopes_bottom_int", "lmm_slopes_top_int")
# GWAS submission script
# submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/BOLT_helpers/1_perform_GWAS_BOLT.sh"
# Filtering submission script
# submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/BOLT_helpers/2_perform_BOLT_filtering.sh"
# Fine-mapping submission script
# submission script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/BOLT_helpers/3_perform_finemapping_BOLT.sh"

submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/BOLT_helpers/3_perform_finemapping_BOLT.sh"

for (s in STRATA_NAMES) {
  for (pr in PARAMETERS) {
    job_options <- paste(
      "-v",
      paste0(
        "STRATA=\"", s, "\",",
        "PARAMETER=\"", pr, "\""
      )
    )
    job_submission <- paste("qsub", job_options, submission_script)
    system(job_submission)
    print(job_submission)
  }
}
