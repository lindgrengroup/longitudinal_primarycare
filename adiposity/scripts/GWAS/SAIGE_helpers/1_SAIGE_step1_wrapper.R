# Author: Samvida S. Venkatesh
# Date: 08/02/22
# Adapted from: Duncan Palmer (https://github.com/astheeggeggs/SAIGE_gene_munging)

# MAKE SURE TO SUBMIT WITH 12 THREADS (or change nThreads in this script)

# library(SAIGE, lib.loc='/well/lindgren/flassen/software/tmp/') 
library(SAIGE)
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--strata", required=TRUE,
                    help = "Which strata we are assessing")
parser$add_argument("--cluster", required=TRUE,
                    help = "Which cluster we are assessing")
parser$add_argument("--adjBaseline",
                    default = "False",
                    help = "Should GWAS be adjusted for baseline trait value?")
parser$add_argument("--plinkFile", 
                    default = "/well/lindgren/UKBIOBANK/ferreira/IMPUTED_association_analysis/grm_files/ukb_cal_v2_qced_pruned")
parser$add_argument("--covars",
                    default = "UKB_assmt_centre,genotyping.array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21",
                    help = "Collection of covariates in the phenotype file to use")
parser$add_argument("--sampleIDColinphenoFile",
                    default = "eid",
                    help = "Sample ID column in the phenotype file")
parser$add_argument("--traitType",
                    default = "binary",
                    help = "What type of trait is the phenotype? Binary or quantitative?")
parser$add_argument("--outputdir",
                    default = "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/SAIGE_results/",
                    help = "Results output directory")

args <- parser$parse_args()
print(args)

# Create some arguments to GLMM from combinations of others
phenoFile <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/traits_for_GWAS/cluster_membership_", 
                    args$strata, ".txt")
outputPrefix <- paste0(args$outputdir, args$strata, "/", args$cluster)
invNormalize <- args$traitType == "quantitative"
covarColList <- unlist(strsplit(args$covars, split=","))
if (grepl("sex_comb", args$strata)) covarColList <- c(covarColList, "sex")
if (args$adjBaseline %in% c("T", "TRUE", "True")) 
  covarColList <- c(covarColList, "baseline_trait")

# Step 1: fit null model
fitNULLGLMM(
  plinkFile = args$plinkFile,
  phenoFile = phenoFile,
  phenoCol = args$cluster,
  covarColList = covarColList,
  sampleIDColinphenoFile = args$sampleIDColinphenoFile,
  traitType = args$traitType,
  invNormalize = invNormalize,
  nThreads = 12, 
  LOCO = FALSE, 
  outputPrefix = outputPrefix,
  IsOverwriteVarianceRatioFile = TRUE)