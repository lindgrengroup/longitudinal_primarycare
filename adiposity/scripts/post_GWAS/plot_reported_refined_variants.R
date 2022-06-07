# Author: Samvida S. Venkatesh
# Date: 20/04/22

library(tidyverse)
theme_set(theme_bw())

# Read in data for cleaning ----

GWAS_zip <- gzfile("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/BOLT_results/BMI_sex_comb_lmm_intercepts_final.txt.gz", "rt")  
GWAS_res <- read.table(GWAS_zip, sep = "\t", header = T, 
                       comment.char = "@", stringsAsFactors = F)

reported_variants <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/BMI_sex_comb/classify_lmm_intercept_variants/reported_snp_list.txt",
                                sep = "\t", header = F,
                                stringsAsFactors = F)$V1

refined_variants <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/BMI_sex_comb/classify_lmm_intercept_variants/refined_variants.txt",
                               sep = "\t", header = F,
                               stringsAsFactors = F)$V1
refined_variants <- gsub(".cma.cojo", "", refined_variants)

# Wrangle data ----

# Prepare columns for cleaning
to_numeric <- c("CHR", "POS", "AF_Tested", "BETA", "SE", "PVALUE")
GWAS_res <- GWAS_res %>% as_tibble() %>%
  mutate(across(all_of(to_numeric), as.numeric)) %>%
  # Flag variants as refined or reported
  mutate(maf = ifelse(AF_Tested < 0.5, AF_Tested, 1-AF_Tested),
         flag = ifelse(SNP %in% reported_variants, "reported",
                       ifelse(SNP %in% refined_variants, "refined", 
                              ifelse(PVALUE < 5E-8, "significant", "ns"))),
         flag = factor(as.character(flag)))

colPalette <- c("#1E3F66", "#FFA500", "#D4D4D4")
names(colPalette) <- c("reported", "refined", "significant")

# Manhattan plots ----

# Format data to plot
datplot <- GWAS_res %>% 
  # get chromosome length
  group_by(CHR) %>% summarise(chr_len = max(POS)) %>%
  # get chromosome position
  mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>% select(-chr_len) %>%
  # add to original results dataset
  left_join(GWAS_res, ., by = c("CHR" = "CHR")) %>%
  # add cumulative position of each SNP
  arrange(CHR, POS) %>% mutate(POS_bp = POS + tot) 

# Axis should just show chromosome number
axisdf <- datplot %>% group_by(CHR) %>% 
  summarise(centre = (max(POS_bp) + min(POS_bp)) / 2)

# Plot
man_BOLT <- ggplot(datplot, aes(x = POS_bp, y = -log10(PVALUE))) +
  geom_point(data = subset(datplot, flag != "ns"), 
             aes(fill = flag), shape = 21, size = 2) +
  scale_fill_manual(values = colPalette, guide = F) +
  geom_point(data = subset(datplot, flag == "ns"),
             aes(color = as.factor(CHR)), shape = 21, alpha = 0.8, size = 1) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 ), guide = "none") +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$centre) +
  scale_y_continuous(limits = c(0, NA)) +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

ggsave(filename = "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/plots/BMI_sex_comb_reported_refined.png",
       plot = man_BOLT,
       device = "png", width = 10, height = 5, units = "in")


