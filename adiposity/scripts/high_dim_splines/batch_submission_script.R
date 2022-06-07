# Author: Samvida S. Venkatesh
# Date: 20/05/2022

# R wrapper to submit sample clustering 
# Need to loop over all of the phenotypes and sex strata, as well as L and M values

PHENOTYPES <- "Weight"
SEX_STRATA <- "sex_comb"
NCLUST <- 5
L <- c(2, 5, 10)
M <- c("random", 1, 2, 5, 10)

submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/2_submit_sample_clustering_scheme.sh"

for (lmin in L) {
  for (myrs in M) {
    job_options <- paste(
      "-v",
      paste0(
        "phenotype=\"", PHENOTYPES, "\",",
        "ss=\"", SEX_STRATA, "\",",
        "nclust=\"", NCLUST, "\",",
        "lmin=\"", lmin, "\",",
        "myrs=\"", myrs, "\""
      )
    )
    job_submission <- paste("qsub", job_options, submission_script)
    system(job_submission)
    print(job_submission)
  }
}
