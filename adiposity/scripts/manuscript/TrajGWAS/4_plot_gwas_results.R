# Author: Samvida S. Venkatesh
# Date: 15/02/22

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

# colour palette: rose, yellow, teal
custom_three_diverge <- c("#D35C79","#C7B241", "#009593")

# Wrangle data for cleaning ----

GWAS_zip <- gzfile(args$inputFile, "rt")  

GWAS_res <- read.csv(GWAS_zip, sep = "\t", header = T, 
                     stringsAsFactors = F)

# Prepare columns for cleaning
to_numeric <- c("chr", "pos", "maf", "infoscore", "hwepval",
                "betapval", "taupval", "jointpval")
GWAS_res <- GWAS_res %>% as_tibble() %>%
  select(c("snpid", "chr", "pos", "varid", "maf", 
           "infoscore", "hwepval",
           "betapval", "betadir", "taupval", "taudir", "jointpval")) %>%
  mutate(across(all_of(to_numeric), as.numeric))

# Cleaning functions ----

# Sanity check: Minor allele frequency
maf_filter <- function (qc_log, dat) {
  res <- dat %>% filter(maf >= 0.01 & maf <= 0.5)
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# SNPs removed with low or implausible MAF < 0.01 or > 1: ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

# Sanity check: INFO score
info_filter <- function (qc_log, dat) {
  res <- dat %>% filter(infoscore >= 0.8)
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# SNPs removed with INFO < 0.8: ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

# Sanity check: Multi-allelic markers
biallelic_filter <- function (qc_log, dat) {
  res <- dat %>% filter(!grepl(";", snpid))
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# multi-allelic SNPs removed: ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

# Remove markers with duplicate entries
duplicate_snps <- function (qc_log, dat) {
  res <- dat %>% 
    distinct(varid, .keep_all = T)
  
  sink(qc_log, append = T)
  cat(paste0("\t", 
             "# SNPs removed due to duplicates: ", 
             nrow(dat) - nrow(res), "\n"))
  sink()
  return (res)
}

# Apply cleaning functions ----

log_file <- args$logFile

cleaned <- maf_filter(log_file, GWAS_res)
cleaned <- info_filter(log_file, cleaned)
cleaned <- biallelic_filter(log_file, cleaned)
cleaned <- duplicate_snps(log_file, cleaned)

# Report metrics post QC
sink(log_file, append = T)
cat(paste0("# SNPs post-QC: ", nrow(cleaned), "\n"))
sink()

# Print file formatted for FUMA and LocusZoom ----

# Sort results by chromosome and position and rename SNP to "chrN:pos"
cleaned <- cleaned %>% arrange(chr, pos) %>%
  mutate(SNP = ifelse(grepl("^rs", snpid), snpid, varid))

to_print <- cleaned[, c("SNP", "chr", "pos", "maf",
                        "betapval", "betadir", "taupval", "taudir", "jointpval")]

# Make sure integers as printed as integers and not in scientific
options(scipen = 999)
write.table(to_print, 
            args$outputFile,
            sep = "\t", row.names = F, quote = F)
options(scipen = 0)

# QQ plots and lambdaGC in each MAF bin ----

maf_max_for_bins <- c(0.01, 0.05, 0.1, 0.51)
NBINS <- length(maf_max_for_bins) - 1

gwas_dat <- cleaned

makeQQ <- function (mafbin, obs_pvals) {
  
  obs_pvals <- obs_pvals[!is.na(obs_pvals) & obs_pvals > 0]
  obs_pvals <- sort(obs_pvals)
  
  # Get genomic inflation factor
  chisq_obs <- qchisq(obs_pvals, df = 1, lower.tail = F)
  lambdaGC <- median(chisq_obs) / qchisq(0.5, 1)
  
  # Get -log10 for plots
  neglog_obs <- -log10(obs_pvals)
  exp_pvals <- (1:length(obs_pvals) - 0.5)/length(obs_pvals)
  neglog_exp <- -log10(exp_pvals)
  
  plot_qq <- data.frame(obs = neglog_obs, exp = neglog_exp)
  
  qq_BOLT <- ggplot(plot_qq, aes(x = neglog_exp, y = neglog_obs)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    labs(title = paste0("MAF bin = ",
                        maf_max_for_bins[mafbin], " - ", 
                        maf_max_for_bins[mafbin+1], 
                        ", lambda = ", lambdaGC),
         x = "Expected -log10(P)", y = "Observed -log10(P)")
  return (qq_BOLT)
}

qq_plots <- lapply(1:NBINS, function (i) {
  sub_gwas <- subset(gwas_dat, gwas_dat$maf >= maf_max_for_bins[i] &
                       gwas_dat$maf < maf_max_for_bins[i+1])
  # Apply to betapvals
  beta_qq <- makeQQ(mafbin = i, obs_pvals = sub_gwas$betapval)
  ggsave(paste0(args$outPlotDir, "beta_pval_qq_mafbin_", i, ".png"),
         beta_qq)
  # Apply to taupvals
  tau_qq <- makeQQ(mafbin = i, obs_pvals = sub_gwas$taupval)
  ggsave(paste0(args$outPlotDir, "tau_pval_qq_mafbin_", i, ".png"),
         tau_qq)
  # Apply to joint pvals
  joint_qq <- makeQQ(mafbin = i, obs_pvals = sub_gwas$jointpval)
  ggsave(paste0(args$outPlotDir, "joint_pval_qq_mafbin_", i, ".png"),
         joint_qq)
})

# Manhattan plots ----

sub_gwas <- gwas_dat

sub_gwas <- sub_gwas %>% 
  # get chromosome length
  group_by(chr) %>% summarise(chr_len = max(pos)) %>%
  # get chromosome position
  mutate(tot = cumsum(chr_len) - chr_len) %>% select(-chr_len) %>%
  # add to original results dataset
  left_join(sub_gwas, ., by = c("chr" = "chr")) %>%
  # add cumulative position of each SNP
  arrange(chr, pos) %>% mutate(BP_pos = pos + tot) %>%
  # plot pval to be min of all pvals
  mutate(plot_pval = pmin(betapval, taupval, jointpval))

# Axis should just show chromosome number
axisdf <- sub_gwas %>% group_by(chr) %>% 
  summarise(centre = (max(BP_pos) + min(BP_pos)) / 2)

# Plot
man_all <- ggplot(sub_gwas, 
                   aes(x = BP_pos, y = -log10(plot_pval))) +
  geom_point(aes(color = as.factor(chr)), alpha = 0.8, size = 1.3) +
  geom_point(data = subset(sub_gwas, -log10(betapval) > 8),
             aes(y = -log10(betapval)),
             color = "#D35C79", shape = 19, size = 2) +
  geom_point(data = subset(sub_gwas, -log10(taupval) > 8),
             aes(y = -log10(taupval)),
             color = "#009593", shape = 19, size = 2) +
  geom_point(data = subset(sub_gwas, -log10(jointpval) > 8),
             aes(y = -log10(jointpval)),
             color = "#C7B241", shape = 19, size = 2) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  scale_x_continuous(label = axisdf$chr, breaks = axisdf$centre) +
  scale_y_continuous(limits = c(0, NA)) +
  theme(legend.position = "none", 
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

ggsave(paste0(args$outPlotDir, "trajgwas_manhattan.png"),
       width = 10, height = 5, units = "in", man_all)

man_beta <- ggplot(sub_gwas, 
                  aes(x = BP_pos, y = -log10(betapval))) +
  geom_point(aes(color = as.factor(chr)), alpha = 0.8, size = 1.3) +
  geom_point(data = subset(sub_gwas, -log10(betapval) > 8),
             aes(y = -log10(betapval)),
             color = "#D35C79", shape = 19, size = 2) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  scale_x_continuous(label = axisdf$chr, breaks = axisdf$centre) +
  scale_y_continuous(limits = c(0, NA)) +
  theme(legend.position = "none", 
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

ggsave(paste0(args$outPlotDir, "beta_manhattan.png"),
       width = 10, height = 5, units = "in", man_beta)

man_tau <- ggplot(sub_gwas, 
                   aes(x = BP_pos, y = -log10(taupval))) +
  geom_point(aes(color = as.factor(chr)), alpha = 0.8, size = 1.3) +
  geom_point(data = subset(sub_gwas, -log10(taupval) > 8),
             aes(y = -log10(taupval)),
             color = "#009593", shape = 19, size = 2) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  scale_x_continuous(label = axisdf$chr, breaks = axisdf$centre) +
  scale_y_continuous(limits = c(0, NA)) +
  theme(legend.position = "none", 
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

ggsave(paste0(args$outPlotDir, "tau_manhattan.png"),
       width = 10, height = 5, units = "in", man_tau)

man_joint <- ggplot(sub_gwas, 
                  aes(x = BP_pos, y = -log10(jointpval))) +
  geom_point(aes(color = as.factor(chr)), alpha = 0.8, size = 1.3) +
  geom_point(data = subset(sub_gwas, -log10(jointpval) > 8),
             aes(y = -log10(jointpval)),
             color = "#C7B241", shape = 19, size = 2) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  scale_x_continuous(label = axisdf$chr, breaks = axisdf$centre) +
  scale_y_continuous(limits = c(0, NA)) +
  theme(legend.position = "none", 
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

ggsave(paste0(args$outPlotDir, "joint_manhattan.png"),
       width = 10, height = 5, units = "in", man_joint)

