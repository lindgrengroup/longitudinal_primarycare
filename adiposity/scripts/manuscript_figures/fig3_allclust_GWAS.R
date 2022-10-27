# Author: Samvida S. Venkatesh
# Date: 22/02/22

library(matrixStats)
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

STRATA <- c("BMI_F", "BMI_M", "BMI_sex_comb",
            "Weight_F", "Weight_M", "Weight_sex_comb")
custom_four_diverge <- c("#D35C79", "#D9AB90", "#9FBCA4", "#009593")
names(custom_four_diverge) <- c("k1", "k2", "k3", "k4")

# Read GWAS results, thinning SNPs as we go ----

# Function to only keep 10,000 SNPs below a threshold to make plotting easier
thinSNPS <- function (pvals, threshold) {
  pass_threshold <- which(pvals < threshold)
  fail_threshold <- which(pvals >= threshold)
  n_thin <- min(length(fail_threshold), 10000)
  keep_failed <- sample(fail_threshold, n_thin, replace = F)
  indices_return <- c(pass_threshold, keep_failed)
  return (indices_return)
}

KEEP_PTHRESH <- 1E-03

gwas_results_thinned <- lapply(STRATA, function (st) {
  res <- lapply(1:4, function (KI) {
    df <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/BOLT_results/",
                            st, "_k", KI, "_final.txt"),
                     sep = "\t", header = T, stringsAsFactors = F, 
                     comment.char = "@")
    res <- df %>%
      distinct(SNP, CHR, POS, PVALUE)
    keepsnps <- thinSNPS(res$PVALUE, KEEP_PTHRESH)
    res <- res[keepsnps, ]
    colnames(res)[which(colnames(res) == "PVALUE")] <- paste0("PVALUE_k", KI,
                                                              "_", st)
    return (res)
  })
  res <- res %>% 
    reduce(full_join, by = c("SNP", "CHR", "POS")) %>%
    distinct()
  return (res)
}) 
dat <- gwas_results_thinned %>% 
  reduce(full_join, by = c("SNP", "CHR", "POS")) %>%
  distinct()

BMI_F_k1 <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/BOLT_results/BMI_F_k1_final.txt",
                       sep = "\t", header = T, stringsAsFactors = F, 
                       comment.char = "@")
BMI_F_k1 <- BMI_F_k1 %>%
  distinct(SNP, CHR, POS, PVALUE)
colnames(BMI_F_k1)[which(colnames(BMI_F_k1) == "PVALUE")] <- "PVALUE_k1_BMI_F"

dat <- full_join(dat, BMI_F_k1, 
                 by = c("SNP", "CHR", "POS", "PVALUE_k1_BMI_F")) %>%
  distinct()

# Plot ----

# Format data to plot
datplot <- dat %>% 
  # get chromosome length
  group_by(CHR) %>% summarise(chr_len = max(POS)) %>%
  # get chromosome position
  mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>% select(-chr_len) %>%
  # add to original results dataset
  left_join(dat, ., by = c("CHR" = "CHR")) %>%
  # add cumulative position of each SNP
  arrange(CHR, POS) %>% mutate(POS_bp = POS + tot)

# Add highlight and annotation information
plot_pvals <- matrixStats::rowMins(as.matrix(datplot %>% select(starts_with("PVALUE"))), 
                                   na.rm = T)
datplot$plot_pval <- plot_pvals 
datplot$highlight <- ifelse(datplot$plot_pval <= 5E-08, "yes", "no")

# Axis should just show chromosome number
axisdf <- datplot %>% group_by(CHR) %>% 
  summarise(centre = (max(POS_bp) + min(POS_bp)) / 2)

# Plot
man_BOLT <- ggplot(datplot, aes(x = POS_bp, y = -log10(plot_pval))) +
  geom_point(data = datplot %>% filter(if_any(starts_with("PVALUE_k1"), ~ . <= 5E-08)), 
             fill = "#D35C79", colour = "#D35C79",
             shape = 19, size = 2) +
  geom_point(data = datplot %>% filter(if_any(starts_with("PVALUE_k2"), ~ . <= 5E-08)), 
             fill = "#D9AB90", colour = "#D9AB90",
             shape = 19, size = 2) +
  geom_point(data = datplot %>% filter(if_any(starts_with("PVALUE_k3"), ~ . <= 5E-08)), 
             fill = "#9FBCA4", colour = "#9FBCA4",
             shape = 19, size = 2) +
  geom_point(data = datplot %>% filter(if_any(starts_with("PVALUE_k4"), ~ . <= 5E-08)), 
             fill = "#009593", colour = "#009593",
             shape = 19, size = 2) +
  geom_point(data = subset(datplot, highlight == "no"),
             aes(color = as.factor(CHR)), shape = 21, alpha = 0.8, size = 1) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22), guide = "none") +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$centre) +
  scale_y_continuous(limits = c(0, NA)) +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

ggsave(filename = paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/plots/all_strata_softclust_manhattan.png"),
       plot = man_BOLT,
       device = "png", width = 10, height = 5, units = "in")
