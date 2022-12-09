# Author: Samvida S. Venkatesh
# Date: 08/03/2022

library(tidyverse)

# R wrapper to submit GP trajectories plotting jobs
# Input:
# 1. path to ids text file: (input format: text file with UKB eid and classification)
# 2. comma-separated list of biomarkers to plot ("all" to run BMI, weight, waist circumference, and WHR)
# 3. name of log file
# 4. directory to where output plots should be stored

GROUPS_OF_INTEREST <- c("sex", "year_of_birth", "baseline_age", 
                        "FUyrs", "FU_n", 
                        "smoking_status", "^ICD_")
CAT_GROUPS <- c("sex", "birth_cohort", "baseline_age_group",
                "fuyrs_group", "fu_n_group",
                "smoking_status")

# Read data ----

PHENOTYPES <- c("BMI", "Weight", "WC", "WHR")

general_covars <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/220131_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, comment.char = "$",
                             stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)
trait_covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[PHENOTYPES]

# Create groups from quantitative data ----

getBirthCohort <- function (x) {
  cut_points <- c(1900, 1940, 1950, 1960, 1971, 2022)
  cut_labels <- c("pre-1940", "'40-'49", "'50-'59", "'60-'70", "post-1970")
  res <- cut(x, breaks = cut_points, labels = cut_labels,
             include.lowest = T, right = F)
  # Ensure the ordering is correct
  res <- factor(as.character(res), levels = cut_labels)
  return (res)
}

getBaselineAgeGroups <- function (x) {
  cut_points <- seq(20, 80, by = 10)
  res <- cut(x, breaks = cut_points, include.lowest = T, right = F)
  return (res)
}

getFUGroups <- function (x) {
  cut_points <- c(0, 4, 10, 15, 1E6)
  cut_labels <- c("0-4", "5-10", "11-15", "16+")
  res <- cut(x, breaks = cut_points, labels = cut_labels,
             include.lowest = T, right = T)
  # Ensure the ordering is correct
  res <- factor(as.character(res), levels = cut_labels)
  return (res)
}

covars <- lapply(trait_covars, function (df) {
  to_calc <- left_join(df, general_covars)
  # Get eid and covariates we want to plot
  to_calc <- to_calc[, c(1, grep(paste0(GROUPS_OF_INTEREST, collapse = "|"), 
                                 colnames(to_calc)))]
  
  to_calc <- to_calc %>% 
    mutate(birth_cohort = getBirthCohort(year_of_birth),
           baseline_age_group = getBaselineAgeGroups(baseline_age),
           fuyrs_group = getFUGroups(FUyrs),
           fu_n_group = getFUGroups(FU_n))
  
  res <- to_calc %>% select(matches(paste0(c("eid", CAT_GROUPS, "^ICD_"),
                                           collapse = "|")))
  return (res)
})

# Loop through phenotypes and categorical groups to submit plotting scripts ----

submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/submit_plot_modelled_trajectories.sh"
loop_groups <- colnames(covars[[1]])[-1]

for (p in PHENOTYPES) {
  for (cat_var in loop_groups) {
    # Write temporary id file
    tmp_file <- covars[[p]][, c("eid", cat_var)]
    write.table(tmp_file, paste0(p, "/tmp_ids_", cat_var, ".txt"),
                sep = "\t", row.names = F, quote = F)
    # Submit job
    job_options <- paste(
      "-v",
      paste0(
        "idFile=\"", paste0(p, "/tmp_ids_", cat_var, ".txt"), "\",",
        "phenotype=\"", p, "\",",
        "outPrefix=\"", paste0(p, "/", cat_var), "\""
      )
    )
    job_submission <- paste("qsub", job_options, submission_script)
    system(job_submission)
    print(job_submission)
  }
}

