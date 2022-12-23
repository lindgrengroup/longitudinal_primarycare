# Author: Samvida S. Venkatesh
# Date: 20/05/2022

# R wrapper to submit pheno and strata-specific jobs
# Need to loop over all of the phenotypes and sex strata

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/scripts/1_submit_hidim_splines.sh"

for (pheno in PHENOTYPES) {
  rvar <- ifelse(pheno == "BMI", 0.5, 5)
  for (ss in SEX_STRATA) {
    job_options <- paste(
      "-v",
      paste0(
        "phenotype=\"", pheno, "\",",
        "sex_strata=\"", ss, "\",",
        "resid_var=\"", rvar, "\""
      )
    )
    job_submission <- paste("qsub", job_options, submission_script)
    system(job_submission)
    print(job_submission)
  }
}
