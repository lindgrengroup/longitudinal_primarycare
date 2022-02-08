# Author: Samvida S. Venkatesh
# Date: 08/02/22
# Adapted from: Duncan Palmer (https://github.com/astheeggeggs/SAIGE_gene_munging)

library(SAIGE, lib.loc='/well/lindgren/flassen/software/tmp/') 
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--strata", required=TRUE,
                    help = "Which strata we are assessing")
parser$add_argument("--cluster", required=TRUE,
                    help = "Which cluster we are assessing")
parser$add_argument("--chr", required = TRUE,
                    help = "Chromosome to run")
parser$add_argument("--bgenDir",
                    default = "/well/ukbb-wtchg/v3/imputation/",
                    help = "Path to directory with bgen and bgi files")
parser$add_argument("--sampleDir",
                    default = "/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/sample_qc/",
                    help = "Path to directory with list of IDs that passed sample QC")
parser$add_argument("--step1ResultsDir", 
                    default = "/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/SAIGE_results/",
                    help = "Path to directory where step 1 results were stored")
parser$add_argument("--step2ResultsDir", 
                    default = "/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/SAIGE_results/",
                    help = "Path to directory where step 2 results should be stored")
args <- parser$parse_args()

# Creat inputs to the SPA test function
bgenFile <- paste0(args$bgenDir, "ukb_imp_chr", args$chr, "_v3.bgen")
bgenFileIndex <- paste0(bgenFile, ".bgi")
sampleFile <- paste0(args$sampleDir, args$strata, "_ids.incl")
GMMATmodelFile <- paste0(args$step1ResultsDir, paste0(args$strata), "/",
                         args$cluster, ".rda")
varianceRatioFile <- paste0(args$step1ResultsDir, paste0(args$strata), "/",
                            args$cluster, ".varianceRatio.txt")
SAIGEOutputFile <- paste0(args$step2ResultsDir, paste0(args$strata), "/",
                          args$cluster, "_chr", args$chr, "_gwas_results.txt")

# Run single-variant analysis

SPAGMMATtest(
  bgenFile = bgenFile,
  bgenFileIndex = bgenFileIndex,
  IsDropMissingDosages = FALSE,
  minMAC = 20,
  sampleFile = sampleFile,
  GMMATmodelFile = GMMATmodelFile,
  varianceRatioFile = varianceRatioFile,
  SAIGEOutputFile = SAIGEOutputFile,
  IsOutputAFinCaseCtrl = TRUE,
  IsOutputNinCaseCtrl = TRUE
  )