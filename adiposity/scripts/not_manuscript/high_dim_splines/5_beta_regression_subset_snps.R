# Author: Samvida S. Venkatesh
# Date: 24/08/22

library(argparse)
library(tidyverse)
library(bigsnpr)
library(betareg)
library(lmtest)
library(broom)

# Parse arguments ----

parser <- ArgumentParser()
parser$add_argument("--phenotype", required=TRUE,
                    help = "Phenotype being tested")
parser$add_argument("--sex_strata", required=TRUE,
                    help = "Sex-strata being tested")
parser$add_argument("--chromosome", required=TRUE,
                    help = "Chromosome being tested")
args <- parser$parse_args()

PHENO <- args$phenotype
SEX_STRATA <- args$sex_strata
CHR <- paste0("chr", args$chromosome)
NCLUSTS <- 4

resdir <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/beta_regression_subset/",
                 PHENO, "_", SEX_STRATA, "/")
dir.create(resdir)

# Read files ----

# Genetic data
geno_fname <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/GWAS/obesity_genotypes/bed_format/",
                       CHR, "_pruned.bed")
tmpfile <- tempfile()
snp_readBed(geno_fname, backingfile = tmpfile)
geno_dat <- snp_attach(paste0(tmpfile, ".rds"))

# Soft clustering probabilities
clust_res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/",
                               PHENO, "_", SEX_STRATA, 
                               "/soft_clustering_probs_", 
                               PHENO, "_", SEX_STRATA, ".txt"),
                        sep = "\t", header = T, stringsAsFactors = F)
clust_res$eid <- as.character(clust_res$eid)

# Covariates
genetic_covars <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/sample_qc/", 
                         PHENO, "_", SEX_STRATA, "_ids_passed_qc.txt"),
                  sep = "\t", header = T)
genetic_covars$eid <- as.character(genetic_covars$IID)

trait_specific_covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[[PHENO]]
trait_specific_covars$eid <- as.character(trait_specific_covars$eid)

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

QCOVARS <- c("baseline_age", "age_sq", 
             "FUyrs", "FU_n",
             paste0("PC", 1:21))
CATCOVARS <- c("sex", "UKB_assmt_centre", 
               "genotyping.array")

if (SEX_STRATA != "sex_comb")
  CATCOVARS <- CATCOVARS[-which(CATCOVARS == "sex")]

# Transform [0,1] probabilities to (0,1) ----

# Based on https://psycnet.apa.org/record/2006-03820-004?doi=1
# Smithson & Verkuilen 2006
squeezeProbs <- function (x, nsamples = 100) {
  squeezed_x <- (x*(nsamples - 1) + 0.5)/nsamples
  return (squeezed_x)
}

# NSAMPLES is based on the number of bootstrap samples used to generate
# the cluster probabilities
squeezed_clust_probs <- apply(clust_res[, -1], 2, FUN = function (x) {
  squeezeProbs(x, nsamples = 100)
})
# ADAPT THIS IF MORE CLUSTERS
squeezed_clust_probs <- data.frame(eid = clust_res$eid, 
                       k1 = squeezed_clust_probs[, "k1"], 
                       k2 = squeezed_clust_probs[, "k2"],
                       k3 = squeezed_clust_probs[, "k3"], 
                       k4 = squeezed_clust_probs[, "k4"]) 

# Create full data ----

genetic_covars <- genetic_covars %>% 
  select(any_of(c("eid", QCOVARS, CATCOVARS)))
trait_specific_covars <- trait_specific_covars %>% 
  select(any_of(c("eid", QCOVARS, CATCOVARS)))
general_covars <- general_covars %>% 
  select(any_of(c("eid", QCOVARS, CATCOVARS)))

all_covars <- inner_join(genetic_covars, trait_specific_covars)
all_covars <- inner_join(all_covars, general_covars) %>%
  mutate(across(all_of(CATCOVARS), function (x) factor(x))) %>%
  mutate(across(all_of(QCOVARS), function (x) as.numeric(x)))

dat_to_model <- inner_join(squeezed_clust_probs, all_covars, 
                           by = "eid")

# Create full data by adding in genotype dosages
Gmat <- as.data.frame(geno_dat$genotypes[])
snpnames <- paste0(geno_dat$map$marker.ID, "_",
                   geno_dat$map$allele1, "_",
                   geno_dat$map$allele2)
colnames(Gmat) <- snpnames
Gmat$eid <- as.character(geno_dat$fam$sample.ID)

full_dat <- left_join(dat_to_model, Gmat, by = "eid") 

# Beta regression ----

performBetaReg <- function (clustk, snpr) {
  # Create formula
  beta_formula <- paste0(clustk, " ~ ", 
                         paste0(c(CATCOVARS, QCOVARS), collapse = " + "), " + ",
                         snpr)
  # Get data to model
  df <- full_dat %>% 
    select(all_of(c("eid", clustk, CATCOVARS, QCOVARS, snpr)))
  df <- df[complete.cases(df), ]
  nsize <- nrow(df)
  
  # Perform beta regression
  model_res <- tryCatch({
    mod <- betareg(formula(beta_formula), data = df)
    res <- tidy(coeftest(mod)) %>%
      filter(term == snpr) %>%
      mutate(N = nsize)
  }, error = function (cond) {
    message(paste0(snpr, " failed due to: "))
    message(cond)
    return (NULL)
  })
  
  return (model_res)
}

# Apply to all clusters ----

lapply(paste0("k", 1:NCLUSTS), function (clustk) {
  cat(paste0("Running cluster: ", clustk, "\n"))
  per_snp <- lapply(snpnames, function (snpr) {
    cat(paste0("\t", "SNP: ", snpr, "\n"))
    res <- performBetaReg(clustk, snpr)
  })
  res <- bind_rows(per_snp) %>%
    rename(SNP = term, beta = estimate, se = std.error,
           zstat = statistic, pval = p.value) %>%
    arrange(pval)
  write.table(res, paste0(resdir, clustk, "_", CHR, ".txt"),
              sep = "\t", row.names = F, quote = F)
})
