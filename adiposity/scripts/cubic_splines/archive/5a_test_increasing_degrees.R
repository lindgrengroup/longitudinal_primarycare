# Author: Samvida S. Venkatesh
# Date: 05/10/21

library(lme4)
library(tidyverse)
library(pheatmap)
theme_set(theme_bw())

set.seed(051021)

# Parse in phenotype argument
args <- commandArgs(trailingOnly = T)
PHENO <- args[1]

SEX_STRATA <- c("F", "M", "sex_comb")

# Read files ----

dat <- 
  readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/TRAINING_SET_adiposity_split_sex.rds")[[PHENO]]
SEX_STRATA <- names(dat)

# Mixed effects model formulas ----

PCs <- paste0("PC", 1:21)
# There is only one data provider for WHR, so this has to be removed from the 
# model
if (PHENO != "WHR") {
  COVARS <- c("data_provider")
} else {
  COVARS <- c()
}

# Function for polynomial spline effect of age
get_formula <- function (n_df_fe, n_df_re, sx) {
  # Adjust for baseline covariates
  if (sx == "sex_comb") { mod_covars <- c(COVARS, "sex") } 
  else { mod_covars <- COVARS }
  
  covars_form <- paste0("value ~ ", 
                        paste0(mod_covars, collapse = " + "), " + ",
                        paste0(PCs, collapse = " + "))
  
  # Add spline effects to formula
  fe_spline <- paste0("poly(age_event, ", n_df_fe, ")")
  re_spline <- paste0("(1 + poly(age_event, ", n_df_re, ") | eid)")
  
  # Full formula
  full_form <- paste0(covars_form, " + ", 
                      fe_spline, " + ", 
                      re_spline)
  return (full_form)
}

# Apply formula to each strata ----

MAX_N_DF <- 10

coefs_matrix <- lapply(SEX_STRATA, function (sx) {
  
  # Data 
  sub_dat <- dat[[sx]]
  
  # Apply formula to 1:MAX_N degrees of freedom, but while there are errors from
  # the random effects structure being too complex, reduce the RE terms; this 
  # will choose the model with most RE terms for the give df
  n_df_mods <- lapply(1:MAX_N_DF, function (n_fe) {
    
    n_re <- n_fe
    mod_formula <- formula(get_formula(n_df_fe = n_fe, n_df_re = n_re, sx))
    res <- try(lmer(mod_formula, sub_dat))
    
    while (class(res) == "try-error" & n_re > 1) {
      n_re <- n_re - 1
      mod_formula <- formula(get_formula(n_df_fe = n_fe, n_df_re = n_re, sx))
      res <- try(lmer(mod_formula, sub_dat))
    }
    
    # Get fixed effect coefficients for spline terms
    spline_coefs <- fixef(res)[grep("poly", names(fixef(res)))]
    names(spline_coefs) <- paste0("poly", 1:n_fe)
    return (spline_coefs)
  })
  coefs_from_mods <- bind_rows(n_df_mods)
  return (coefs_from_mods)
})
names(coefs_matrix) <- SEX_STRATA

saveRDS(coefs_matrix, 
        paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/spline_coefficients_", 
               PHENO, ".rds"))

# Plot coefficients matrix ----

pdf(paste0("plots/increasing_spline_degrees_", PHENO, ".pdf"), onefile = T)
plot_coefs_matrix <- lapply(SEX_STRATA, function (sx) {
  dat_plot <- abs(coefs_matrix[[sx]])
  pheatmap(dat_plot, 
           cluster_rows = F, cluster_cols = F,
           scale = "none", treeheight_col = 0,
           color = colorRampPalette(c("white", "red"))(50),
           labels_row = paste0("n_degrees = ", 1:MAX_N_DF),
           labels_col = paste0("term", 1:MAX_N_DF),
           main = paste(PHENO, sx), fontsize = 8)
})
dev.off()
