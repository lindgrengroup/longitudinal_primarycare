get_control_parameters <- function() {
  control <- list()
  
  ##########################################
  # Data and working directories
  ##########################################
  if (Sys.info()["user"] == "qgj547") { #George
    control$adiposity_root <- "/well/lindgren-ukbb/projects/ukbb-11867/georgenicholson/github_repos/longitudinal_primarycare/adiposity"
    # control$biobank_completion_root <- "/well/lindgren-ukbb/projects/ukbb-11867/georgenicholson/github_repos/longitudinal_primarycare/biobank_completion"
  } else { # Samvida
    control$adiposity_root <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity"
    # control$biobank_completion_root <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/biobank_completion"
  }
  # setwd(control$biobank_completion_root)
  dir.create("plots", showWarnings = FALSE)
  dir.create("output", showWarnings = FALSE)
  
  ##########################################
  # Some global parameters
  ##########################################
  control$PHENO <- c("BMI", "Weight", "WHR", "WC")[1]
  control$SEX_STRATA <- c("F", "M", "sex_comb")[3]
  control$do_plots <- FALSE
  
  ##########################################
  # Age binning
  ##########################################
  control$AGE_BIN_WID <- 1
  control$MIN_AGE <- 21
  control$MAX_AGE <- 79
  
  ##########################################
  # Trait binning
  ##########################################
  control$N_Y_BIN <- 100
  
  ##########################################
  # Define spline basis 
  ##########################################
  control$AR1_RHO <- 0.975 # Main smoothness parameter, typically useful between about .9 and .995
  control$AR1_NOISE_SD <- 5 # Can leave this as is typically
  control$AR1_INTERCEPT_SD <- 100 # A big number for noninformative intercept in AR1
  control$RESID_NOISE_VAR <- .01
  control$NDF_SPLINE <- 200 # DF of spline
  
  ##############################################################################
  # Parameters for learning local autocorrelation structure 
  ##############################################################################
  control$N_GRAD_CLASSES <- 5 # Must be an odd number
  control$N_MC_SAMPLES <- 1000
  control$NTASKS_01 <- 1000
  
  ##############################################################################
  # Parameters for HMM
  ##############################################################################
  control$HMM_RESID_NOISE_VAR <- .05
  control$N_LIK_EVAL_PTS <- 20
  control$N_ITS <- 1000
  control$NTASKS_02 <- 1000
  
  ##############################################################################
  # Parameters for optimising variance parameters 
  ##############################################################################
  control$N_SUBJ_PER_OPTIM <- 200
  control$NTASKS_03 <- 100

  ##############################################################################
  # Parameters for fitting genetic assoc model
  ##############################################################################
  control$NTASKS_04 <- 1000
  control$NTASKS_04a <- 50
  
  return(control)
  
}