# Author: Samvida S. Venkatesh
# Date: 30/06/21

library(tidyverse)
theme_set(theme_bw())

# Read files ----

dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/visually_QCd_adiposity.rds")
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/model_covariates.rds")

PHENO <- names(dat)
PCs <- paste0("PC", 1:21)
COVARS <- c("baseline_age", "age_sq", "FUyrs", "FU_n", "PCs", 
            "data_provider")
SEX_STRATA <- c("F", "M", "sex_comb")

# Subset to white ancestry for these analyses ----

dat <- lapply(PHENO, function (ph) {
  res <- merge(dat[[ph]], covars[[ph]], by = "eid")
  res <- subset(res, res$ancestry == "white")
  return (res)
})
names(dat) <- PHENO

# Calculate variance explained by covariates ----

## Function to calculate R2 per covariate ----

var_explained_per_covar <- function (ph, sx) {
  
  # Subset the right data and covariates
  if (sx == "sex_comb") {
    MOD_COVARS <- c(COVARS, "sex")
    sub_dat <- dat[[ph]]
  } else {
    MOD_COVARS <- COVARS
    sub_dat <- subset(dat[[ph]], dat[[ph]]$sex == sx)
  }
  # Calculate R2 per covariate (all genetic PCs together)
  per_covar_R2 <- lapply(MOD_COVARS, function (cov) {
    if (cov == "PCs") {
      mod_form <- paste0("value ~ ", paste0(PCs, collapse = " + "))
    } else {
      mod_form <- paste0("value ~ ", cov)
    }
    mod_R2 <- summary(lm(as.formula(mod_form), 
                         data = sub_dat))$r.squared
    return (mod_R2)
  })
  names(per_covar_R2) <- MOD_COVARS
  per_covar_R2 <- bind_cols(per_covar_R2)
  per_covar_R2$trait <- ph
  per_covar_R2$sex_strata <- sx
  
  return (per_covar_R2)
} 

## Run function and plot results per sex strata -----

all_res <- lapply(SEX_STRATA, function (sx) {
  all_pheno <- lapply(PHENO, function (ph) {
    return (var_explained_per_covar(ph, sx))
  })
  # Return results table
  all_pheno <- bind_rows(all_pheno)
  # Return results plot
  plot_df <- pivot_longer(all_pheno, cols = c(-trait, -sex_strata),
                          names_to = "covariate",
                          values_to = "R2")
  covar_plot <- ggplot(plot_df, aes(x = trait, y = R2,
                                    fill = covariate)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(palette = "Set3") +
    labs(title = paste0("Proportion of trait variance explained in sex: ", 
                        sx)) +
    coord_flip() 
  return (list(res_tbl = all_pheno,
          res_plot = covar_plot))
})

## Write table of results ----

all_tbls <- bind_rows(lapply(all_res, function (x) x[["res_tbl"]] ))
write.table(all_tbls, "var_explained/per_covar_R2.txt",
            sep = "\t", quote = F, row.names = F)

## Print plots ----

pdf("var_explained/per_covar_R2.pdf", onefile = T)
print(lapply(all_res, function (x) x[["res_plot"]] ))
dev.off()

# Full covariate models ----

## Run full model and add residuals and fitted value to data -----

adj_dat <- lapply(PHENO, function (ph) {
  # Log file for variance explained
  log_var <- paste0("var_explained/all_covar_R2_", ph, ".txt")
  
  per_sex <- lapply(SEX_STRATA, function (sx) {
    # Subset the right data and covariates
    if (sx == "sex_comb") {
      MOD_COVARS <- c(COVARS, "sex")
      sub_dat <- dat[[ph]]
    } else {
      MOD_COVARS <- COVARS
      sub_dat <- subset(dat[[ph]], dat[[ph]]$sex == sx)
    }
    # Define model
    mod_form <- paste0("value ~ ", 
                       paste0(MOD_COVARS[which(MOD_COVARS != "PCs")], 
                              collapse = " + "), " + ",
                       paste0(PCs, collapse = " + "))
    mod <- lm(as.formula(mod_form), data = sub_dat)
    
    # Create new dataframe with residuals and fitted values
    res_df <- sub_dat
    res_df$adj_value <- mod$residuals
    res_df$fitted_covs <- fitted(mod)
    
    # Log variance explained
    sink(log_var, append = T)
    cat(paste0("Proportion of ", 
               ph, " variance explained by all covariates in ", sx),
        "\n")
    print(summary(mod)$r.squared)
    cat("\n")
    sink()
    
    return (res_df)
  })
  names(per_sex) <- SEX_STRATA
  # Save adjusted values within each trait
  saveRDS(per_sex, paste0("adj_traits/adj_", ph, ".rds"))
  return (per_sex)
})

