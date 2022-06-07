# Author: Samvida S. Venkatesh
# Date: 01/03/2022

# R wrapper to submit qctool jobs
# Loop over all lead variants in each parameter (across all strata)

submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/4_extract_dosages.sh"

# PARAMETERS <- c("lmm_intercepts", "lmm_slopes_adj_baseline", "cspline_intercepts")
PARAMETERS <- "Weight_sex_comb_clusters"
lead_variants <- lapply(PARAMETERS, function (pr) {
  df <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/lead_snps/",
                          pr, "_ALL_variants.txt"),
                   sep = "\t", header = T, stringsAsFactors = F)
  return (df)
})
names(lead_variants) <- PARAMETERS

for (pr in PARAMETERS) {
  df <- lead_variants[[pr]]
  for (vi in 1:nrow(df)) {
    job_options <- paste(
      "-v",
      paste0(
        "VARID=\"", df$rsid[vi], "\",",
        "CHR=\"", df$chromosome[vi], "\""
      )
    )
    job_submission <- paste("qsub", job_options, submission_script)
    system(job_submission)
    print(job_submission)
  }
}
