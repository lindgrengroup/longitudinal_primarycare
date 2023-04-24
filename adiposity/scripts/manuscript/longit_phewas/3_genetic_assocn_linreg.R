# Author: Samvida S. Venkatesh
# Date: 11/10/21

library(tidyverse)

# Read data ----

infile_path <- "" # REDACTED
gpdat_path <- "" # REDACTED
gen_resources_path <- "" # REDACTED
outfile_path <- "" # REDACTED

PHENOTYPES <- read.table(paste0(infile_path, "/code_lists/qcd_traits_available.txt"),
                         sep = "\t", header = F, stringsAsFactors = F)$V1
REMOVE_PHENOS <- c("BMI", "Weight", "WC", "WHR", "FAI", "Progesterone", "Prolactin")
PHENOTYPES <- PHENOTYPES[!PHENOTYPES %in% REMOVE_PHENOS]
SEX_STRATA <- c("F", "M", "sex_comb")

blups <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SEX_STRATA, function (sx) {
    df <- read.table(paste0(infile_path, "/longit_phewas/lmm_models/",
                            p, "_", sx, "_blups_full_model.txt"), 
                     sep = "\t", header = T, stringsAsFactors = F)
    df$eid <- as.character(df$eid)
    colnames(df) <- c("eid", "b0", "b1")
    return (df)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(blups) <- PHENOTYPES

# IDs that passed genotyping QC
ids_passed_geno_qc <- lapply(PHENOTYPES, function (p) {
  df <- read.table(paste0(infile_path, "/longit_phewas/sample_qc/", 
                          p, "_ids_passed_qc.txt"), sep = "\t", header = T, stringsAsFactors = F)
  df$eid <- as.character(df$eid)
  return (df)
})
names(ids_passed_geno_qc) <- PHENOTYPES

# IDs in discovery GWAS
gp_ids <- read.table(paste0(gpdat_path, "/data/all_ids_in_discovery_gp_ukb_dat.txt"),
                     sep = "\t", header = F, stringsAsFactors = F)$V1
gp_ids <- as.character(gp_ids)

# Genotypes / dosages at lead SNPs
vars_to_replicate <- read.table(paste0(gpdat_path, "/data/lead_snps_to_replicate.txt"),
                                sep = "\t", header = T, stringsAsFactors = F)
VARIDS <- vars_to_replicate$SNP
# For variants without rsids, build variant name
rename_vars <- grep("^chr", VARIDS)
VARIDS[rename_vars] <- paste0(vars_to_replicate$SNP[rename_vars], "_",
                              vars_to_replicate$Tested_Allele[rename_vars], "_", 
                              vars_to_replicate$Other_Allele[rename_vars])
VARIDS <- gsub("^chr", "", VARIDS)

# Genotypes / dosages at variants of interest
var_dosages <- lapply(VARIDS, function (varid) {
  res <- read.table(paste0(gpdat_path, "/sample_variant_counts/",
                           varid, "_dosages.txt"),
                    sep = " ", header = T, stringsAsFactors = F)
  # Remove first row, which contains info on type of column and columns 
  # 2, 3, 4 (ID repeat, missingness, sex)
  res <- res[-1, c(1, 5)]
  colnames(res) <- c("eid", varid)
  return (res)
})
names(var_dosages) <- VARIDS

covars <- readRDS(paste0(gpdat_path, "/data/covariates.rds"))[PHENOTYPES]

general_covars <- read.table(paste0(gen_resources_path, "/220504_QCd_demographic_covariates.txt"),
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

ADJ_COVARS <- c("baseline_age", "age_sq", "year_of_birth",
                "FU_n", "FUyrs", "UKB_assmt_centre", "genotyping.array",
                paste0("PC", 1:21)) # add sex for sc analyses

# Convert rs429358 dosages to genotypes (hard call) ----

# Convert dosages to 0/1/2 genotypes based on threshold
DOSAGE_THRESHOLD_0 <- 0.5
DOSAGE_THRESHOLD_2 <- 1.5
var_dosages_hardcall <- lapply(VARIDS, function (varid) {
  res <- var_dosages[[varid]] %>% 
    mutate(genotype = ifelse(!!as.symbol(varid) < DOSAGE_THRESHOLD_0, "0",
                             ifelse(!!as.symbol(varid) > DOSAGE_THRESHOLD_2, 
                                    "2", "1")),
           eid = as.character(eid)) %>%
    dplyr::select(all_of(c("eid", "genotype")))
  return (res)
})
names(var_dosages_hardcall) <- VARIDS

# Run regressions ----

getEffectSizePretty <- function (adj_covars, dat) {
  # Only run the regression if there are > 100 samples, otherwise we will likely be underpowered anyway 
  if (nrow(dat) > 100) {
    # RINT b1
    trait_to_rint <- dat$b1
    dat$rinted_b1 <- qnorm((rank(trait_to_rint) - 0.5) / sum(!is.na(trait_to_rint)))
    
    # Get formula for modelling SNP effect
    model_formula <- formula(paste0("rinted_b1 ~ genotype + b0 + ", 
                                    paste0(adj_covars, collapse = " + ")))
    mod_res <- lm(formula = model_formula, data = dat)
    
    # Return effect size, S.E., P-value (of SNP)
    to_return <- as.data.frame(t(summary(mod_res)$coefficients["genotype", c(1,2,4)]))
    colnames(to_return) <- c("beta", "se", "pvalue")
  } else {
    to_return <- NULL
  }
  return (to_return)
}

all_res <- lapply(VARIDS, function (varid) {
  cat(paste0("Running SNP: ", varid, "\n")) 
  
  per_var <- lapply(PHENOTYPES, function (p) {
    per_sex <- lapply(SEX_STRATA, function (sx) {
      print(paste0(p, "_", sx))
      
      lmm_df <- blups[[p]][[sx]]
      # Only retain the IDs that passed genotyping sample QC
      lmm_df <- lmm_df %>% filter(eid %in% ids_passed_geno_qc[[p]]$eid)
      
      # Add covariates
      lmm_df <- left_join(lmm_df, covars[[p]], by = "eid") 
      lmm_df <- left_join(lmm_df, general_covars, by = "eid")
      lmm_df <- left_join(lmm_df, ids_passed_geno_qc[[p]], by = "eid")
      
      # Add genotype
      lmm_df <- left_join(lmm_df, var_dosages_hardcall[[varid]], by = "eid") %>%
        mutate(genotype = as.numeric(genotype))
      
      to_adj <- ADJ_COVARS
      if (sx == "sex_comb") to_adj <- c(ADJ_COVARS, "sex")
      
      # Get effect size for full dataset
      full_dat_res <- getEffectSizePretty(adj_covars = to_adj,
                                          dat = lmm_df)
      if (!is.null(full_dat_res)) {
        full_dat_res <- full_dat_res %>%
          mutate(n = nrow(lmm_df),
                 discovery_ids_incl = T)
      }
      
      # Get effect size within group that does not have any weight or BMI discovery data
      keep_ids <- which(!lmm_df$eid %in% gp_ids)
      sub_dat_res <- getEffectSizePretty(adj_covars = to_adj,
                                         dat = lmm_df[keep_ids, ])
      if (!is.null(sub_dat_res)) {
        sub_dat_res <- sub_dat_res %>%
          mutate(n = length(keep_ids),
                 discovery_ids_incl = F)
      }
      
      to_write <- bind_rows(full_dat_res, sub_dat_res) %>%
        mutate(varid = varid,
               sex_strata = sx,
               phenotype = p)
      return (to_write)
    })
    names(per_sex) <- SEX_STRATA
    per_sex <- bind_rows(per_sex)
    return (per_sex)
  })
  names(per_var) <- PHENOTYPES
  per_var <- bind_rows(per_var)
  
  write.table(per_var,
              paste0(outfile_path, "/longit_phewas/", varid, "_all_results.txt"),
              sep = "\t", quote = F, row.names = F)
  return (per_var)
})
all_res <- bind_rows(all_res)

