# Author: Samvida S. Venkatesh
# Date: 25/05/21

library(tidyverse)
theme_set(theme_bw())

# Read files to combine ----

args <- commandArgs(trailingOnly = T)
PHENO <- args[1]
STRATA <- args[2]

SEXES <- c("F", "M", "sexcomb")

GWAS_res <- lapply(SEXES, function (s) {
  res <- readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/GWAS/BOLT_filtered/", 
                        PHENO, "_", s, "_", STRATA, ".rds"))
  return (res)
})
names(GWAS_res) <- SEXES

# Get only the most significant SNPs in each file ----

top_SNPs <- lapply(GWAS_res, function (r) {
  return (r$SNP[r$P_BOLT_LMM_INF < 5e-04])
})
top_SNPs <- unique(unlist(top_SNPs))

GWAS_comb <- lapply(SEXES, function (s) {
  res <- subset(GWAS_res[[s]], GWAS_res[[s]]$SNP %in% top_SNPs)
  res <- res[, c("SNP", "CHR", "BP", "BETA", "SE", "P_BOLT_LMM_INF")]
  colnames(res)[6] <- "PVAL"
  res <- res %>% distinct(SNP, CHR, BP, .keep_all = T)
  res$sex <- s
  return (res)
})
GWAS_comb <- bind_rows(GWAS_comb)
GWAS_comb <- GWAS_comb %>% group_by(SNP, CHR, BP) %>%
  mutate(LOWEST_PVAL = min(PVAL)) %>% ungroup()
GWAS_comb <- pivot_wider(GWAS_comb, id_cols = c(SNP, CHR, BP, LOWEST_PVAL),
                         names_from = sex,
                         values_from = c(BETA, SE, PVAL))

# Plot results ----

dat <- GWAS_comb %>% 
  # get chromosome length
  group_by(CHR) %>% summarise(chr_len = max(BP)) %>%
  # get chromosome position
  mutate(tot = cumsum(chr_len) - chr_len) %>% select(-chr_len) %>%
  # add to original results dataset
  left_join(GWAS_comb, ., by = c("CHR" = "CHR")) %>%
  # add cumulative position of each SNP
  arrange(CHR, BP) %>% mutate(BP_pos = BP + tot) %>%
  # Add highlight and annotation information
  mutate(highlight = ifelse(-log10(PVAL_sexcomb) > 8 | 
                              (-log10(PVAL_F) > 8 & -log10(PVAL_M) > 8), 
                            "sex-comb", 
                            ifelse(-log10(PVAL_F) > 8, "fem_only",
                                   ifelse(-log10(PVAL_M) > 8, 
                                          "male_only", "no"))))

# Axis should just show chromosome number
axisdf <- dat %>% group_by(CHR) %>% 
  summarise(centre = (max(BP_pos) + min(BP_pos)) / 2)

# THINK OF A BETTER WAY TO DO THIS

# Plot
man_BOLT <- ggplot(dat, aes(x = BP_pos, y = -log10(LOWEST_PVAL))) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22)) +
  geom_point(data = subset(dat, highlight != "no"),
             aes(fill = highlight), size = 2) +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$centre) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = STRATA) +
  theme(legend.position = "none", 
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/plots/GWAS/combined_", 
           PHENO, "_", STRATA, ".pdf"), onefile = T)
print(man_BOLT)
dev.off()