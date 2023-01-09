# Author: Samvida S. Venkatesh
# Date: 22/02/22

library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

args <- commandArgs(trailingOnly = T)
STRATA <- args[1]

# Function to thin SNPs as results are read (save memory) ----

# Function to only keep 100,000 SNPs below a threshold to make plotting easier
thinSNPS <- function (pvals, threshold) {
  pass_threshold <- which(pvals < threshold)
  fail_threshold <- which(pvals >= threshold)
  n_thin <- min(length(fail_threshold), 100000)
  keep_failed <- sample(fail_threshold, n_thin, replace = F)
  indices_return <- c(pass_threshold, keep_failed)
  return (indices_return)
}

# Read GWAS results, thinning SNPs as we go ----

KEEP_PTHRESH <- 1E-03
GWAS_zip <- gzfile(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/SAIGE_results/baseline_age_adj/",
                          STRATA, "/", STRATA, "_k1_final.txt.gz"), "rt")  
k1_full <- read.table(GWAS_zip,
                  sep = "\t", header = T, stringsAsFactors = F, 
                  comment.char = "@")
k1_full$cluster <- "k1"

gwas_results_thinned <- lapply(2:5, function (KI) {
  GWAS_zip <- gzfile(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/SAIGE_results/baseline_age_adj/",
                            STRATA, "/", STRATA, "_k", KI, "_final.txt.gz"), 
                     "rt")
  df <- read.table(GWAS_zip,
                   sep = "\t", header = T, stringsAsFactors = F, 
                   comment.char = "@")
  # Only keep SNPs with pval < KEEP_PTHRESH 
  # Thin SNPs
  res <- df %>% filter(PVALUE < KEEP_PTHRESH)
  # Add cluster identity
  res$cluster <- paste0("k", KI)
  return (res)
})
dat <- bind_rows(gwas_results_thinned)
dat <- bind_rows(k1_full, dat)

# Plot ----

colPalette <- brewer.pal(5, "Set1")
names(colPalette) <- paste0("k", 1:5)

dat$cluster <- factor(dat$cluster)

# Format data to plot
datplot <- dat %>% 
  # get chromosome length
  group_by(CHR) %>% summarise(chr_len = max(POS)) %>%
  # get chromosome position
  mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>% select(-chr_len) %>%
  # add to original results dataset
  left_join(dat, ., by = c("CHR" = "CHR")) %>%
  # add cumulative position of each SNP
  arrange(CHR, POS) %>% mutate(POS_bp = POS + tot) %>%
  # Add highlight and annotation information
  mutate(highlight = ifelse(PVALUE < 5E-8, "yes", "no"))

# Axis should just show chromosome number
axisdf <- datplot %>% group_by(CHR) %>% 
  summarise(centre = (max(POS_bp) + min(POS_bp)) / 2)

# Plot
man_BOLT <- ggplot(datplot, aes(x = POS_bp, y = -log10(PVALUE))) +
  geom_point(aes(color = as.factor(CHR)), shape = 21, alpha = 0.8, size = 1) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 ), guide = "none") +
  geom_point(data = subset(datplot, highlight == "yes"), 
             aes(fill = cluster), shape = 21, size = 2) +
  scale_fill_manual(values = colPalette, guide = F) +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$centre) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(title = STRATA) +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

ggsave(filename = paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/plots/baseline_age_adj/",
                         STRATA, "/all_clusters_", STRATA, ".png"),
       plot = man_BOLT,
       device = "png", width = 10, height = 5, units = "in")
