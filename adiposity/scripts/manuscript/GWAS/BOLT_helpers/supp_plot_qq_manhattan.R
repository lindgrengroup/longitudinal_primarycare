# Author: Samvida S. Venkatesh
# Date: 21/05/21

library(argparse)
library(tidyverse)
theme_set(theme_bw())

# Read in arguments ----

parser <- ArgumentParser()
parser$add_argument("--strata", required=TRUE,
                    help = "Strata")
parser$add_argument("--parameter", required=TRUE,
                    help = "Cluster (combination)")
args <- parser$parse_args()

main_filepath <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/GWAS/")

STRATA <- args$strata
PARAMETER <- args$parameter

# Wrangle data ----

gwas_dat <- read.table(paste0(main_filepath, "BOLT_results/", STRATA, "_", PARAMETER, "_final.txt"), 
                       sep = "\t", header = T, 
                       comment.char = "@", stringsAsFactors = F)

# Prepare columns for cleaning
to_numeric <- c("CHR", "POS", "AF_Tested", "BETA", "SE", "PVALUE")
gwas_dat <- gwas_dat %>% as_tibble() %>%
  mutate(across(all_of(to_numeric), as.numeric)) %>%
  mutate(maf = ifelse(AF_Tested < 0.5, AF_Tested, 1-AF_Tested))

# QQ plots and lambdaGC in each MAF bin ----

maf_max_for_bins <- c(0.01, 0.05, 0.1, 0.51)
NBINS <- length(maf_max_for_bins) - 1

qq_plots <- lapply(1:NBINS, function (i) {
  sub_gwas <- subset(gwas_dat, gwas_dat$maf >= maf_max_for_bins[i] &
                       gwas_dat$maf < maf_max_for_bins[i+1])
  # Get genomic inflation factor
  lambdaGC <- median(qchisq(1 - sub_gwas$PVALUE, 1)) / qchisq(0.5, 1)
  
  obs_pvals <- sub_gwas$PVALUE
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
  ggsave(paste0(main_filepath, "plots/", STRATA, "/", PARAMETER, "/qq_mafbin_", i, ".png"),
         qq_BOLT)
})

# Manhattan plots ----

# # Function to only keep 10,000 SNPs below a threshold to make plotting easier
# thinSNPS <- function (pvals, threshold) {
#   pass_threshold <- which(pvals < threshold)
#   fail_threshold <- which(pvals >= threshold)
#   n_thin <- min(length(fail_threshold), 10000)
#   keep_failed <- sample(fail_threshold, n_thin, replace = F)
#   indices_return <- c(pass_threshold, keep_failed)
#   return (indices_return)
# }
# 
# # keep only some SNPs to reduce plot size
# sub_gwas <- gwas_dat[thinSNPS(gwas_dat$P_BOLT_LMM_INF, 1E-03), ]

sub_gwas <- gwas_dat

sub_gwas <- sub_gwas %>% 
  # get chromosome length
  group_by(CHR) %>% summarise(chr_len = max(POS)) %>%
  # get chromosome position
  mutate(tot = cumsum(chr_len) - chr_len) %>% select(-chr_len) %>%
  # add to original results dataset
  left_join(sub_gwas, ., by = c("CHR" = "CHR")) %>%
  # add cumulative position of each SNP
  arrange(CHR, POS) %>% mutate(BP_pos = POS + tot) %>%
  # Add highlight and annotation information
  mutate(highlight = ifelse(-log10(PVALUE) > 8, "yes", "no"))

# Axis should just show chromosome number
axisdf <- sub_gwas %>% group_by(CHR) %>% 
  summarise(centre = (max(BP_pos) + min(BP_pos)) / 2)

# Plot
man_BOLT <- ggplot(sub_gwas, aes(x = BP_pos, y = -log10(PVALUE))) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
  geom_point(data = subset(sub_gwas, highlight == "yes"), 
             color = "orange", size = 2) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$centre) +
  scale_y_continuous(limits = c(0, NA)) +
  theme(legend.position = "none", 
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

ggsave(paste0(main_filepath, "plots/", STRATA, "/", PARAMETER, "/manhattan.png"),
       width = 10, height = 5, units = "in", man_BOLT)

