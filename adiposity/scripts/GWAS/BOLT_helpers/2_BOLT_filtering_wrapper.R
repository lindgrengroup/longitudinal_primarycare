# Author: Samvida S. Venkatesh
# Date: 21/05/21

library(argparse)
library(tidyverse)
theme_set(theme_bw())

# Read in arguments ----

parser <- ArgumentParser()
parser$add_argument("--inputFile", required=TRUE,
                    help = "Path to combined and filtered summary statistics")
parser$add_argument("--logFile", required=TRUE,
                    help = "Path to log file to store # SNPs cleaned")
parser$add_argument("--outputFile", required = TRUE,
                    help = "Output file name for filtered summary statistics")
parser$add_argument("--outPlotDir", required = TRUE,
                    help = "Path to store output Manhattan and QQ-plots")
args <- parser$parse_args()

# Wrangle data for cleaning ----

GWAS_zip <- gzfile(args$inputFile, "rt")  
GWAS_res <- read.table(GWAS_zip, sep = "\t", header = T, 
                       comment.char = "@", stringsAsFactors = F)

# Prepare columns for cleaning
to_numeric <- c("CHR", "BP", "F_MISS", "A1FREQ", "BETA", "SE", "F_MISS",
                "CHISQ_BOLT_LMM_INF", "P_BOLT_LMM_INF", "P_LINREG")
GWAS_res <- GWAS_res %>% as_tibble() %>%
  mutate(across(all_of(to_numeric), as.numeric)) %>%
  mutate(maf = ifelse(A1FREQ < 0.5, A1FREQ, 1-A1FREQ))

# Cleaning functions ----

# Sanity check: minor allele frequency
maf_clean <- function (qc_log, dat) {
  res <- dat %>% filter(maf > 0 & maf <= 0.5)
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# SNPs removed with implausible MAF (<0 or >1): ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

# Sanity check: multi-allelic markers
biallelic_filter <- function (qc_log, dat) {
  res <- dat %>% filter(!grepl(";", SNP))
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# multi-allelic SNPs removed: ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

# Missingness < 5%
fmiss <- function (qc_log, dat) {
  res <- dat %>% filter(F_MISS < 0.05)
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# SNPs removed with missingness > 5%: ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

# Implausibly large standard error (> 10)
extreme_effect <- function (qc_log, dat) {
  res <- dat %>% filter(SE < 10)
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# SNPs removed with extreme standard error > 10: ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

# Remove markers with duplicate entries
duplicate_snps <- function (qc_log, dat) {
  res <- dat %>% 
    mutate(dup_check = paste0(CHR, ":", BP, "_", ALLELE1, "_", ALLELE0)) %>%
    distinct(dup_check, .keep_all = T)
  
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# SNPs removed due to duplicates: ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

# Apply cleaning functions ----

log_file <- args$logFile

cleaned <- maf_clean(log_file, GWAS_res)
cleaned <- biallelic_filter(log_file, cleaned)
cleaned <- fmiss(log_file, cleaned)
cleaned <- extreme_effect(log_file, cleaned)
cleaned <- duplicate_snps(log_file, cleaned)

# Report metrics post QC
sink(log_file, append = T)
cat(paste0("\t", "# SNPs post-QC: ", nrow(cleaned), "\n"))
sink()

# Print file formatted for FUMA and LocusZoom ----

# Sort results by chromosome and position and rename SNP to "chrN:pos"
cleaned <- cleaned %>% arrange(CHR, BP) %>%
  mutate(SNP = ifelse(grepl("^rs", SNP), SNP, 
                      paste0("chr", sub("_.*", "", SNP))))

to_print <- cleaned[, c("SNP", "CHR", "BP", 
                        "ALLELE1", "ALLELE0", "A1FREQ", 
                        "BETA", "SE", "P_BOLT_LMM_INF")]
colnames(to_print) <- c("SNP", "CHR", "POS", 
                        "Tested_Allele", "Reference_Allele", "AF_Tested",
                        "BETA", "SE", "PVALUE")
# Make sure integers as printed as integers and not in scientific
options(scipen = 999)
write.table(to_print, 
            args$outputFile,
            sep = "\t", row.names = F, quote = F)
options(scipen = 0)

# QQ plots and lambdaGC in each MAF bin ----

gwas_dat <- cleaned

maf_max_for_bins <- c(0.01, 0.05, 0.1, 0.51)
NBINS <- length(maf_max_for_bins) - 1

qq_plots <- lapply(1:NBINS, function (i) {
  sub_gwas <- subset(gwas_dat, gwas_dat$maf >= maf_max_for_bins[i] &
                       gwas_dat$maf < maf_max_for_bins[i+1])
  # Get genomic inflation factor
  lambdaGC <- median(sub_gwas$CHISQ_BOLT_LMM_INF) / qchisq(0.5, 1)
  
  obs_pvals <- sub_gwas$P_BOLT_LMM_INF
  obs_pvals <- obs_pvals[!is.na(obs_pvals) & obs_pvals > 0]
  obs_pvals <- sort(obs_pvals)
  neglog_obs <- -log10(obs_pvals)
  exp_pvals <- (1:length(obs_pvals) - 0.5)/length(obs_pvals)
  neglog_exp <- -log10(exp_pvals)
  
  plot_qq <- data.frame(obs = neglog_obs, exp = neglog_exp)
  
  qq_BOLT <- ggplot(plot_qq, aes(x = neglog_exp, y = neglog_obs)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    labs(title = paste0("MAF bin = ",
                        maf_max_for_bins[i], " - ", 
                        maf_max_for_bins[i+1], 
                        ", lambda = ", lambdaGC),
         x = "Expected -log10(P)", y = "Observed -log10(P)")
  ggsave(paste0(args$outPlotDir, "qq_mafbin_", i, ".png"),
         qq_BOLT)
})

# Manhattan plots ----

# Function to only keep 10,000 SNPs below a threshold to make plotting easier
thinSNPS <- function (pvals, threshold) {
  pass_threshold <- which(pvals < threshold)
  fail_threshold <- which(pvals >= threshold)
  n_thin <- max(length(fail_threshold), 10000)
  keep_failed <- sample(fail_threshold, n_thin, replace = F)
  indices_return <- c(pass_threshold, keep_failed)
  return (indices_return)
}

# keep only some SNPs to reduce plot size
sub_gwas <- gwas_dat[thinSNPS(gwas_dat$P_BOLT_LMM_INF, 1E-03), ]

sub_gwas <- sub_gwas %>% 
  # get chromosome length
  group_by(CHR) %>% summarise(chr_len = max(BP)) %>%
  # get chromosome position
  mutate(tot = cumsum(chr_len) - chr_len) %>% select(-chr_len) %>%
  # add to original results dataset
  left_join(sub_gwas, ., by = c("CHR" = "CHR")) %>%
  # add cumulative position of each SNP
  arrange(CHR, BP) %>% mutate(BP_pos = BP + tot) %>%
  # Add highlight and annotation information
  mutate(highlight = ifelse(-log10(P_BOLT_LMM_INF) > 8, "yes", "no"))

# Axis should just show chromosome number
axisdf <- sub_gwas %>% group_by(CHR) %>% 
  summarise(centre = (max(BP_pos) + min(BP_pos)) / 2)

# Plot
man_BOLT <- ggplot(sub_gwas, aes(x = BP_pos, y = -log10(P_BOLT_LMM_INF))) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
  geom_point(data = subset(sub_gwas, highlight == "yes"), 
             color = "orange", size = 2) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$centre) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(title = STRATA) +
  theme(legend.position = "none", 
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

ggsave(paste0(args$outPlotDir, "manhattan.png"),
       width = 10, height = 5, units = "in", man_BOLT)

