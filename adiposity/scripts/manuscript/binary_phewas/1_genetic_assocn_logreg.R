# Author: Samvida S. Venkatesh
# Date: 11/10/21

library(tidyverse)

# Read data ----

gpdat_path <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity" # REDACTED
gen_resources_path <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources" # REDACTED
outfile_path <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/binary_phewas" # REDACTED

# ID-phenotype matrix
eid_pheno_matrix <- read.table(paste0(gen_resources_path, "/eid_phenotype_matrix.txt"), 
                               sep = "\t", header = T, stringsAsFactors = F)
rownames(eid_pheno_matrix) <- eid_pheno_matrix$eid
eid_pheno_matrix <- eid_pheno_matrix[, -1]

# Disease dictionary
dictionary <- read.table(paste0(gen_resources_path, "/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt"),
                         sep = "\t", header = T, stringsAsFactors = F,
                         comment.char = "@", quote = "")
DISEASES <- dictionary$phenotype
colnames(eid_pheno_matrix) <- DISEASES

# Genotypes / dosages at rs429358
apoe_dosages <- read.table(paste0(gpdat_path, "/sample_variant_counts/rs429358_dosages.txt"),
                           sep = " ", header = T, stringsAsFactors = F)
# Remove first row, which contains info on type of column and columns 
# 2, 3, 4 (ID repeat, missingness, sex)
apoe_dosages <- apoe_dosages[-1, c(1, 5)]
colnames(apoe_dosages) <- c("eid", "dosage")

general_covars <- read.table(paste0(gen_resources_path, "/220504_QCd_demographic_covariates.txt"),
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

SEX_STRATA <- c("F", "M", "sex_comb")

# Convert rs429358 dosages to genotypes (hard call) ----

# Convert dosages to 0/1/2 genotypes based on threshold
DOSAGE_THRESHOLD_0 <- 0.5
DOSAGE_THRESHOLD_2 <- 1.5
apoe_geno <- apoe_dosages %>% 
  mutate(rs429358 = ifelse(dosage < DOSAGE_THRESHOLD_0, 0,
                           ifelse(dosage > DOSAGE_THRESHOLD_2, 
                                  2, 1)),
         eid = as.character(eid)) %>%
  dplyr::select(all_of(c("eid", "rs429358")))

eid_pheno_matrix$eid <- as.character(rownames(eid_pheno_matrix))
full_dat <- inner_join(apoe_geno, eid_pheno_matrix, by = "eid")
full_dat$sex <- general_covars$sex[match(full_dat$eid, general_covars$eid)]

# Run regressions ----

getEffectSizePretty <- function (dat, sx = "sex_comb") {
  # Only run the regression if there are > 100 cases, otherwise we will likely be underpowered anyway 
  dat <- dat[complete.cases(dat), ]
  colnames(dat)[which(colnames(dat) %in% DISEASES)] <- "disease"
  if (sum(dat$disease) > 100) {
    # Get formula for modelling SNP effect
    if (sx == "sex_comb") {
      model_formula <- formula("disease ~ rs429358 + sex")
    } else {
      dat <- dat %>% filter(sex == sx)
      model_formula <- formula("disease ~ rs429358")
    }

    mod_res <- glm(formula = model_formula, family = binomial(link = "logit"),
                   data = dat)
    
    # Return effect size, S.E., P-value (of SNP)
    to_return <- as.data.frame(t(summary(mod_res)$coefficients["rs429358", c(1,2,4)]))
    colnames(to_return) <- c("beta", "se", "pvalue")
  } else {
    to_return <- NULL
  }
  return (to_return)
}

all_res <- lapply(DISEASES, function (dx) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    print(paste0(dx, "_", sx))
    
    # Subset to data
    dat_model <- full_dat %>%
      select(all_of(c("eid", "sex", "rs429358", dx)))
    
    to_write <- NULL
    if (nrow(dat_model) > 0) {
      # Get effect size
      full_dat_res <- getEffectSizePretty(dat = dat_model, sx = sx)
      if (!is.null(full_dat_res)) {
        to_write <- full_dat_res %>%
          mutate(sex_strata = sx,
                 disease = dx)
      }
    } 
    return (to_write)
  })
  names(per_sex) <- SEX_STRATA
  per_sex <- bind_rows(per_sex) 
  return (per_sex)
})
names(all_res) <- DISEASES
all_res <- bind_rows(all_res)

# convert to OR and LCI/UCI
all_res <- all_res %>%
  mutate(OR = exp(beta),
         lci = exp(beta - 1.96*se),
         uci = exp(beta + 1.96*se))

write.table(all_res,
            paste0(outfile_path, "/rs429358_all_results.txt"),
            sep = "\t", quote = F, row.names = F)

