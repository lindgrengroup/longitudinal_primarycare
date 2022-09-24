# Author: Samvida S. Venkatesh
# Date: 22/09/22

library(GRAB) 
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--strata", required=TRUE,
                    help = "Which strata we are assessing")
parser$add_argument("--chr", required = TRUE,
                    help = "Chromosome to run")
parser$add_argument("--bgenDir",
                    default = "/well/ukbb-wtchg/v3/imputation/",
                    help = "Path to directory with bgen and bgi files")
parser$add_argument("--sampleFile",
                    default = "/well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample",
                    help = "List of sample IDs in bgen files")
parser$add_argument("--step1ResultsDir", 
                    default = "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/POLMM_results/",
                    help = "Path to directory where step 1 results were stored")
parser$add_argument("--step2ResultsDir", 
                    default = "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/POLMM_results/",
                    help = "Path to directory where step 2 results should be stored")
args <- parser$parse_args()

# Creat inputs to the single-variant test function ----

bgenFile <- paste0(args$bgenDir, "ukb_imp_chr", args$chr, "_v3.bgen")
bgenFileIndex <- paste0(bgenFile, ".bgi")
OutputFile <- paste0(args$step2ResultsDir, paste0(args$strata), "/",
                          "chr", args$chr, "_gwas_results.txt")

objPOLMMFile <- paste0(args$step1ResultsDir, args$strata, 
                       "_step1_results.RData")
load(objPOLMMFile)

# Run single-variant analysis ----

GRAB.Marker(obj.POLMM, 
            GenoFile = bgenFile,
            GenoFileIndex = c(bgenFileIndex, args$sampleFile),
            OutputFile = OutputFile)
