# Author: Samvida S. Venkatesh
# Date: 19/07/22

library(tidyverse)
theme_set(theme_bw())

filepath_main <- "" # REDACTED

# Read data ----

STRATA <- c("BMI_F", "BMI_M", "BMI_sex_comb",
            "Weight_F", "Weight_M", "Weight_sex_comb")

all_lead_snps <- lapply(STRATA, function (s) {
  res <- read.table(paste0(filepath_main,
                           "lead_snps/", s, "_b0_final.lead_snps.txt"),
                    sep = "\t", header = T, stringsAsFactors = F)
  return (res)
})
names(all_lead_snps) <- STRATA

novel_refined_snps <- lapply(STRATA, function (s) {
  res <- read.table(paste0(filepath_main, s, "/to_annotate_results_refined_novel_snps.txt"),
                    sep = "\t", header = T, stringsAsFactors = F)
  return (res)
})
names(novel_refined_snps) <- STRATA

# Get percent of lead SNPs that are novel or refined ----

unique_lead_snps <- lapply(STRATA, function (s) {
  res <- all_lead_snps[[s]] %>% 
    mutate(chrpos = paste0(chromosome, ":", position)) %>%
    select(all_of(c("rsid", "chrpos"))) 
  return (res)
})
unique_lead_snps <- bind_rows(unique_lead_snps) %>%
  distinct()

unique_novel_refined_snps <- lapply(STRATA, function (s) {
  res <- novel_refined_snps[[s]] %>% 
    rename(rsid = SNP_og, chr = CHR_og, pos = POS_og) %>%
    mutate(chrpos = paste0(chr, ":", pos)) %>%
    select(all_of(c("rsid", "chrpos", "status")))
  return (res)
})
unique_novel_refined_snps <- bind_rows(novel_refined_snps) %>%
  distinct()

# Calculate proportion of variance explained by novel or refined SNPs ----

col_palette <- c("#009593", "#C7B241", "#D35C79")
names(col_palette) <- c("published", "refined", "novel_refined")

prop_var_novel_refined <- lapply(STRATA, function (s) {
  df <- all_lead_snps[[s]]
  df <- df %>%
    mutate(varexp = (beta^2)*2*maf*(1 - maf),
           status = ifelse(rsid %in% novel_refined_snps[[s]]$SNP_og, 
                           "novel_refined", "published"))
  vexp <- df %>% group_by(status) %>%
    summarise(varexp = sum(varexp), 
              n = n()) %>%
    mutate(strata = s)
  return (vexp)
})
prop_var_novel_refined <- bind_rows(prop_var_novel_refined) %>%
  mutate(strata = factor(strata, levels = c("Weight_sex_comb", "Weight_M", "Weight_F",
                                            "BMI_sex_comb", "BMI_M", "BMI_F")))

bmi_plot <- ggplot(prop_var_novel_refined[grep("BMI", prop_var_novel_refined$strata), ], 
                   aes(x = strata, y = varexp)) +
  geom_col(aes(fill = status), width = 0.7) +
  scale_fill_manual(values = col_palette, guide = "none") +
  scale_y_continuous(limits = c(0, 0.04)) +
  coord_flip() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

tiff(filename = paste0(filepath_main, "bmi_varexpl.tiff"),
     height = 5, width = 5, units = "cm",
     res = 300)
print(bmi_plot)
dev.off()

weight_plot <- ggplot(prop_var_novel_refined[grep("Weight", prop_var_novel_refined$strata), ], 
                   aes(x = strata, y = varexp)) +
  geom_col(aes(fill = status), width = 0.7) +
  scale_fill_manual(values = col_palette, guide = "none") +
  scale_y_continuous(limits = c(0, 0.075)) +
  coord_flip() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

tiff(filename = paste0(filepath_main, "weight_varexpl.tiff"),
     height = 5, width = 5, units = "cm",
     res = 300)
print(weight_plot)
dev.off()


# Grab sumstats for novel SNPs in each strata ----

novel_snps <- lapply(STRATA, function (s) {
  df <- all_lead_snps[[s]]
  novel_snps <- novel_refined_snps[[s]] %>% filter(status == "novel")
  novel_snps <- novel_snps$SNP_og
  res <- df %>% filter(rsid %in% novel_snps)
  res <- res %>% 
    mutate(strata = s,
           chrpos = paste0(chromosome, ":", position, ":", allele1, ":", allele2),
           beta_se = paste0(signif(beta, 3), " (", signif(se, 3), ")"))
  return (res[, c("strata", "rsid", "chrpos", "maf", "beta_se", "p")])
})
novel_snps <- bind_rows(novel_snps)

