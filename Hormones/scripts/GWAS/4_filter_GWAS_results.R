# Author: Samvida S. Venkatesh
# Date: 21/05/21

# Submit with 2 cores

library(tidyverse)
theme_set(theme_bw())

set.seed(210521)

# Read files ----

SNPs_passed_QC <- read.table("/well/lindgren/UKBIOBANK/samvida/full_primary_care/GWAS/snps_passed_QC_211102/passed_QC.txt",
                             sep = "\t", header = F, stringsAsFactors = F)

args <- commandArgs(trailingOnly = T)
STRATA <- args[1]

filename <- paste0("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS/BOLT_results/",
                   STRATA, "_assoc.stats.gz")

log_file <- paste0("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS/log_files/", 
                   STRATA, "_post_GWAS_filtering.txt")

# Cleaning functions ----

find_maf <- function (dat) {
  # Calculate which of the alleles is the minor allele and flip as necessary
  res <- dat %>% mutate(MINOR_ALLELE = ifelse(A1FREQ < 0.5, ALLELE1, ALLELE0),
                        MAF = ifelse(A1FREQ < 0.5, A1FREQ, 1 - A1FREQ))
  return (res)
}

hwe_computed <- function (qc_log, dat) {
  # Previously computed list of SNPs with 
  # HWE pval > 1E-06, biallelic
  res <- dat %>% filter(SNP %in% SNPs_passed_QC$V1)
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# SNPs removed with HWE pval < 1E-06 and non-biallelic: ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

maf_filter <- function (qc_log, dat) {
  # BOLT should have already applied these filters but 
  # double-check
  res <- dat %>%
    filter(MAF > 0.001)
  sink(qc_log, append = T)
  cat(paste0("\t",
             "# SNPs removed with MAF < 0.1%: ",
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

extreme_effect <- function (qc_log, dat) {
  # SNPs with implausibly large standard error (> 10)
  res <- dat %>% filter(SE < 10)
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# SNPs removed with extreme standard error > 10: ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

duplicate_snps <- function (qc_log, dat) {
  tmp <- dat
  tmp$marker_name <- paste0(tmp$CHROM, ":", tmp$GENPOS, ":", 
                            tmp$ALLELE1, tmp$ALLELE0)
  tmp$dups <- duplicated(tmp$marker_name, fromLast = F) | 
    duplicated(tmp$marker_name, fromLast = T)
  res <- tmp %>% filter(!dups)
  
  res <- res[, colnames(dat)]
  
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# SNPs removed due to duplicates: ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

# Apply cleaning functions ----

if (file.exists(filename)) {
  
  GWAS_zip <- gzfile(filename, "rt")  
  GWAS_res <- read.table(GWAS_zip, sep = "\t", header = T, 
                         comment.char = "~", stringsAsFactors = F)
  
  # Report metrics
  sink(log_file, append = T)
  cat(paste0("** Phenotype and Strata: ", STRATA, "\n",
             "\t", "# SNPs pre-QC: ", nrow(GWAS_res), "\n"))
  sink()
  
  GWAS_res$CHR <- as.numeric(GWAS_res$CHR)
  GWAS_res$BP <- as.numeric(GWAS_res$BP)
  GWAS_res$F_MISS <- as.numeric(GWAS_res$F_MISS)
  GWAS_res$A1FREQ <- as.numeric(GWAS_res$A1FREQ)
  GWAS_res$BETA <- as.numeric(GWAS_res$BETA)
  GWAS_res$SE <- as.numeric(GWAS_res$SE)
  GWAS_res$CHISQ_BOLT_LMM_INF <- as.numeric(GWAS_res$CHISQ_BOLT_LMM_INF)
  GWAS_res$P_BOLT_LMM_INF <- as.numeric(GWAS_res$P_BOLT_LMM_INF)
  GWAS_res$P_LINREG <- as.numeric(GWAS_res$P_LINREG)
  
  cleaned <- hwe_computed(log_file, GWAS_res)
  cleaned <- find_maf(cleaned)
  cleaned <- maf_filter(log_file, cleaned)
  cleaned <- extreme_effect(log_file, cleaned)
  cleaned <- duplicate_snps(log_file, cleaned)
  
  # Report metrics post QC
  sink(log_file, append = T)
  cat(paste0("\t", "# SNPs post-QC: ", nrow(cleaned), "\n"))
  sink()
  
  write.table(cleaned, 
              paste0("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS/BOLT_filtered/", 
                     STRATA, ".txt"),
              sep = "\t", row.names = F, quote = F)
}

# QQ plots and lambdaGC in each MAF bin ----

gwas_dat <- read.table(paste0("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS/BOLT_filtered/", 
                              STRATA, ".txt"),
                       sep = "\t", header = T)

maf_max_for_bins <- c(0.001, 0.01, 0.05, 0.1, 0.5)
NBINS <- length(maf_max_for_bins) - 1

qq_plots <- lapply(1:NBINS, function (i) {
  sub_gwas <- subset(gwas_dat, gwas_dat$MAF > maf_max_for_bins[i] &
                       gwas_dat$MAF <= maf_max_for_bins[i+1])
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
  ggsave(paste0("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS/plots/qq_", 
                STRATA, "mafbin_", i, ".png"),
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

ggsave(paste0("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS/plots/manhattan_", 
              STRATA, ".png"),
       width = 10, height = 5, units = "in", man_BOLT)

