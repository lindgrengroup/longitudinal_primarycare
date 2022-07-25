# Author: Samvida S. Venkatesh
# Date: 19/07/22

library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

# Read data ----

STRATA <- c("BMI_F", "BMI_M", "BMI_sex_comb",
            "Weight_F", "Weight_M", "Weight_sex_comb")

known_gwascat_snps <- read.table("/well/lindgren/samvida/Resources/GWASCatalog/gwascat_obesity_associations_hg19.bed",
                                 header = F, sep = "\t")
colnames(known_gwascat_snps) <- c("CHR", "POS0", "POS1", "SNP")
known_gwascat_snps <- known_gwascat_snps %>%
  mutate(CHRPOS = paste0(CHR, ":", POS1)) %>%
  select(all_of(c("SNP", "CHRPOS")))

# Classification of SNPs as novel or refined
novel_refined_snps <- lapply(STRATA, function (s) {
  res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/", 
                           s, "/classify_lmm_intercept_variants/annotated_results_refined_novel_snps.txt"),
                    sep = "\t", header = T, stringsAsFactors = F)
  res <- res %>% 
    mutate(strata = s) %>%
    rename(SNP = SNP_og, CHRPOS = CHRPOS_og) %>%
    select(all_of(c("SNP", "CHRPOS", "status")))
  return (res)
})
names(novel_refined_snps) <- STRATA

# Only read full GWAS results for weight sex-comb and BMI sex-comb
# For the rest, we can read in just the lead SNP results as these are the ones
# that we will be classifying as reported/refined/novel

bmi_file <- gzfile("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/BOLT_results/BMI_sex_comb_lmm_intercepts_final.txt.gz",
                   "rt")
weight_file <- gzfile("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/BOLT_results/Weight_sex_comb_lmm_intercepts_final.txt.gz",
                      "rt")

background_gwas_bmi <- read.table(bmi_file,
                                  sep = "\t", header = T, stringsAsFactors = F, 
                                  comment.char = "@")
background_gwas_bmi <- background_gwas_bmi %>%
  mutate(strata = "BMI_sex_comb") %>%
  select(all_of(c("SNP", "CHR", "POS", "PVALUE", "strata")))

background_gwas_weight <- read.table(weight_file,
                                     sep = "\t", header = T, stringsAsFactors = F, 
                                     comment.char = "@")
background_gwas_weight <- background_gwas_weight %>%
  mutate(strata = "Weight_sex_comb") %>%
  select(all_of(c("SNP", "CHR", "POS", "PVALUE", "strata")))

background_gwas <- bind_rows(background_gwas_bmi,
                             background_gwas_weight)

# Read in lead SNP results for all GWAS and classify as reported, novel, refined
all_lead_snps <- lapply(STRATA, function (s) {
  res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/lead_snps/", s, "_lmm_intercepts_final.lead_snps.txt"),
                    sep = "\t", header = T, stringsAsFactors = F)
  res <- res %>% 
    mutate(strata = s) %>%
    rename(SNP = rsid, CHR = chromosome, POS = position, PVALUE = p) %>%
    select(all_of(c("SNP", "CHR", "POS", "PVALUE", "strata")))
  
  res$status <- novel_refined_snps[[s]]$status[match(res$SNP, 
                                                     novel_refined_snps[[s]]$SNP)]
  res$status[is.na(res$status)] <- "reported"
  
  return (res)
})
names(all_lead_snps) <- STRATA
all_lead_snps <- bind_rows(all_lead_snps)

# Prep plotting data ----

# For background GWAS, colour SNPs with previous GWAS catalog associations
# as well as SNPs with P < 5E-8

# Colour palette
# light grey, light blue, navy, teal green, amber, rose
col_palette <- c("#EBEBEB", "#D4F4FA", "#005580", 
                 "#009593", "#C7B241", "#D35C79")
names(col_palette) <- c("odd_nonsig", "even_nonsig", 
                        "sig", 
                        "reported", "refined", "novel")

background_gwas <- background_gwas %>%
  mutate(status = ifelse(SNP %in% known_gwascat_snps$SNP & PVALUE < 5e-8, "reported",
                               ifelse(CHR %% 2 == 0 & PVALUE >= 5e-8, "even_nonsig",
                                      ifelse(CHR %% 2 != 0 & PVALUE >= 5e-8, "odd_nonsig",
                                                    ifelse(PVALUE < 5e-8, "sig", NA)))))

plot_dat <- bind_rows(background_gwas, all_lead_snps)
plot_dat$status <- as.factor(plot_dat$status)

# Format data to plot
for_plot <- plot_dat %>% 
  # get chromosome length
  group_by(CHR) %>% summarise(chr_len = max(POS)) %>%
  # get chromosome position
  mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>% select(-chr_len) %>%
  # add to original results dataset
  left_join(plot_dat, ., by = c("CHR" = "CHR")) %>%
  # add cumulative position of each SNP
  arrange(CHR, POS) %>% mutate(POS_bp = POS + tot) 

# Axis should just show chromosome number
axisdf <- for_plot %>% group_by(CHR) %>% 
  summarise(centre = (max(POS_bp) + min(POS_bp)) / 2)

# Plot
manhattan_plot <- ggplot(for_plot, 
                         aes(x = POS_bp, y = -log10(PVALUE)),
                         fill = status, colour = status) +
  geom_point(data = for_plot %>% filter(!status %in% c("novel", "refined", "reported")),
             aes(fill = status, colour = status), shape = 19, size = 1) +
  geom_point(data = for_plot %>% filter(status == "reported"), 
             aes(fill = status, colour = status), shape = 19, size = 1) +
  geom_point(data = for_plot %>% filter(status == "refined"), 
             aes(fill = status, colour = status), shape = 19, size = 1) +
  geom_point(data = for_plot %>% filter(status == "novel"), 
             aes(fill = status, colour = status), shape = 19, size = 1) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
  scale_colour_manual(values = col_palette, guide = "none") +
  scale_fill_manual(values = col_palette, guide = "none") +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$centre) +
  scale_y_continuous(limits = c(0, NA)) +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

ggsave(filename = "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/manhattan_novel_refined_reported_intercept_variants.png",
       plot = manhattan_plot,
       device = "png", width = 10, height = 5, units = "in")

