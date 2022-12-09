# Author: Samvida S. Venkatesh
# Date: 14/09/22

library(argparse)
library(GRAB)

parser <- ArgumentParser()
parser$add_argument("--strata", required=TRUE,
                    help = "Which strata we are assessing")
parser$add_argument("--covars",
                    default = "baseline_age,age_sq,UKB_assmt_centre,genotyping.array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21",
                    help = "Collection of covariates in the phenotype file to use")
parser$add_argument("--sampleIDColinphenoFile",
                    default = "eid",
                    help = "Sample ID column in the phenotype file")
parser$add_argument("--ordinalColinphenoFile",
                    default = "clust",
                    help = "Ordinal phenotype column in the phenotype file")
parser$add_argument("--outputdir",
                    default = "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/POLMM_results/",
                    help = "Results output directory")

args <- parser$parse_args()
print(args)

# Read in data ----

genoFile <- "/well/lindgren-ukbb/projects/ukbb-11867/ferreira/IMPUTED_association_analysis/grm_files/ukb_cal_v2_qced_pruned.bed"
sparseGRMFile <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/GRAB/sparseGRM.txt"

phenoFile <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/traits_for_GWAS/", 
                    args$strata, "_assigned_clusts.txt")
outputPrefix <- paste0(args$outputdir, args$strata)

gwas_dat <- read.table(phenoFile, sep = "\t", header = T, stringsAsFactors = F)

# Ensure that ordinal phenotype column is a factor and ordered correctly
ordPhenoLevels <- sort(unique(gwas_dat[, args$ordinalColinphenoFile]))
gwas_dat$phenotype <- factor(gwas_dat[, args$ordinalColinphenoFile],
                             levels = ordPhenoLevels)

# Create formula ----

covarColList <- unlist(strsplit(args$covars, split=","))
if (grepl("sex_comb", args$strata)) covarColList <- c(covarColList, "sex")

gwas_formula <- formula(paste0("phenotype ~ ", 
                               paste0(covarColList, collapse = " + ")))

# Double-check that all samples are in genotype files  ----

geno_ids <- read.table(paste0(gsub(".bed", "", genoFile), ".fam"),
                       header = F, stringsAsFactors = F)
colnames(geno_ids)[1] <- "IID"

gwas_dat <- gwas_dat[gwas_dat$eid %in% geno_ids$IID, ]
 
# Step 1: fit null model ----

obj.POLMM = GRAB.NullModel(formula = gwas_formula,
                           data = gwas_dat, 
                           subjData = gwas_dat[, args$sampleIDColinphenoFile], 
                           method = "POLMM", 
                           traitType = "ordinal",
                           GenoFile = genoFile,
                           SparseGRMFile = sparseGRMFile,
                           control = list(showInfo = FALSE, 
                                          LOCO = FALSE, 
                                          tolTau = 0.2, 
                                          tolBeta = 0.1))

save(obj.POLMM, file = paste0(outputPrefix, "_step1_results.RData")) 
