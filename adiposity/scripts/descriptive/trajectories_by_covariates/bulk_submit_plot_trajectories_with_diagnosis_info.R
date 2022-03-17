# Author: Samvida S. Venkatesh
# Date: 08/02/2022

# R wrapper to submit GP trajectories plotting jobs
# Input:
# 1. comma-separated list of diagnoses to plot ("all" to run all diseases in /well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt)
# 2. name of biomarker to plot
# 3. name of log file
# 4. directory where output plots should be stored

submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/scripts/submit_plot_trajectories_with_diagnosis_info.sh"
DIAGS_TO_PLOT <- "all"
OUTPLOTDIR <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/testing_plot_traj/diagnosis_trajectories/"

BIOMARKERS <- c("BMI", "Weight")
BM_FU <- c(20, 20)
names(BM_FU) <- BIOMARKERS

for (bm in BIOMARKERS) {
  logfile_name <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/testing_plot_traj/all_diags_",
                         bm, "_log.txt")
  job_options <- paste(
    "-v",
    paste0(
      "diagnoses=\"", DIAGS_TO_PLOT, "\",",
      "biomarker=\"", bm, "\",",
      "followUpYears=\"", BM_FU[[bm]], "\",",
      "logFile=\"", logfile_name, "\",",
      "outPlotDir=\"", OUTPLOTDIR, "\""
    )
  )
  job_submission <- paste("qsub", job_options, submission_script)
  system(job_submission)
  print(job_submission)
}
