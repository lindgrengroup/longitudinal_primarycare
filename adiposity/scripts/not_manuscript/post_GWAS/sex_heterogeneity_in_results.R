# Author: Samvida S. Venkatesh
# Date: 24/10/2022

library(argparse)
library(tidyverse)
theme_set(theme_bw())

# Get arguments ----

parser <- ArgumentParser()
parser$add_argument("--phenotype", required = TRUE,
                    help = "Phenotype")
parser$add_argument("--model_term", required = TRUE,
                    help = "Model term: lmm slope or intercept")

args <- parser$parse_args()

PHENO <- args$phenotype
TERM <- args$model_term

# Read data ----

gzf <- gzfile(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/BOLT_results/maf_1/",
                     PHENO, "_F_", TERM, "_final.txt.gz"), "rt")
gzm <- gzfile(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/BOLT_results/maf_1/",
                     PHENO, "_M_", TERM, "_final.txt.gz"), "rt")

dat_f <- read.table(gzf, sep = "\t", header = T, stringsAsFactors = F)
colnames(dat_f)[6:9] <- c("AF_Tested_F", "BETA_F", "SE_F", "PVALUE_F")

dat_m <- read.table(gzm, sep = "\t", header = T, stringsAsFactors = F)
colnames(dat_m)[6:9] <- c("AF_Tested_M", "BETA_M", "SE_M", "PVALUE_M")

# Combine data across sexes and perform Z-test ----

full_dat <- inner_join(dat_f, dat_m)

full_dat <- full_dat %>% 
  mutate(het_zstat = (BETA_F - BETA_M)/sqrt(SE_F^2 + SE_M^2),
         het_pval = pnorm(het_zstat, 0, 1, lower.tail = T))

# Get summary for results ----

# Plot qq-plot for het-pvals

obs_pvals <- full_dat$het_pval
obs_pvals <- obs_pvals[!is.na(obs_pvals) & obs_pvals > 0]
obs_pvals <- sort(obs_pvals)
neglog_obs <- -log10(obs_pvals)
exp_pvals <- (1:length(obs_pvals) - 0.5)/length(obs_pvals)
neglog_exp <- -log10(exp_pvals)

plot_qq <- data.frame(obs = neglog_obs, exp = neglog_exp)

qqres <- ggplot(plot_qq, aes(x = neglog_exp, y = neglog_obs)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Expected -log10(P)", y = "Observed -log10(P)")
ggsave(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/plots/",
              PHENO, "_", TERM, "_sex_heterogeneity_pvals.png"),
       qqres)

write.table(full_dat[full_dat$het_pval <= 5E-08, ], 
            paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/BOLT_results/",
                   PHENO, "_", TERM, "_sex_heterogeneity_results.txt"),
            sep = "\t", row.names = F, quote = F)


