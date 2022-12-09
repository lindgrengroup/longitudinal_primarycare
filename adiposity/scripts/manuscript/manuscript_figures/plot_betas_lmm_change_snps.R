# Author: Samvida S. Venkatesh
# Date: 08/09/2022

library(tidyverse)
theme_set(theme_bw())

# colour palette: rose, teal, grey
custom_three_diverge <- c("#D35C79","#009593", "#666666")
names(custom_three_diverge) <- c("F", "M", "sex_comb")

# resdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2211_models/GWAS/change_snps/"
# 
# # Read data ----
# 
# PHENOTYPES <- c("BMI", "Weight")
# SEX_STRATA <- c("F", "M", "sex_comb")
# 
# change_snps_sumstats <- lapply(PHENOTYPES, function (p) {
#   res <- lapply(SEX_STRATA, function (sx) {
#     df_intercepts <- read.table(paste0(resdir, p, "_", sx, "_intercepts.txt"),
#                             sep = "\t", header = F, stringsAsFactors = F)
#     colnames(df_intercepts) <- c("SNP", "CHR", "POS", 
#                                 "Tested_Allele", "Other_Allele",
#                                 "AF_Tested", "BETA", "SE", "PVALUE")
#     df_intercepts$term <- "intercept"
#     
#     df_slopes <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/change_snps/",
#                                       p, "_", sx, "_slopes.txt"),
#                                sep = "\t", header = F, stringsAsFactors = F)
#     colnames(df_slopes) <- c("SNP", "CHR", "POS", 
#                                 "Tested_Allele", "Other_Allele",
#                                 "AF_Tested", "BETA", "SE", "PVALUE")
#     df_slopes$term <- "slope"
#     
#     df <- bind_rows(df_intercepts, df_slopes) %>%
#       mutate(phenotype = p,
#              sex_strata = sx,
#              strata = paste0(p, "_", sx))
#     return (df)
#   })
#   res <- bind_rows(res)
#   return (res)
# })
# change_snps_sumstats <- bind_rows(change_snps_sumstats)
# 
# write.table(change_snps_sumstats,
#             paste0(resdir, "all_strata_sumstats.txt"), sep = "\t",
#             quote = F, row.names = F)

# Plot results ----

change_snps_sumstats <- read.table("C:/Users/samvida/Documents/Lindgren Group/Adiposity_Primary_Care/2211_models/rs429358_all_lmm_terms.txt",
                                   sep = "\t", header = T, stringsAsFactors = F)

# Flip effect allele
plot_dat <- change_snps_sumstats %>%
  mutate(v = "rs429358",
         sex_strata = factor(sex_strata, levels = c("F", "M", "sex_comb")),
         term = factor(term, levels = c("b0", "b1")),
         strata = factor(paste0(phenotype, "_", sex_strata), 
                         levels = c("Weight_sex_comb", "Weight_M",
                                    "Weight_F",
                                    "BMI_sex_comb", "BMI_M", "BMI_F")),
         BETA_FLIP = -BETA,
         uci = BETA_FLIP + 1.96*SE,
         lci = BETA_FLIP - 1.96*SE,
         sig_lty = factor(ifelse(PVALUE < 5E-8, "yes", "no"),
                          levels = c("yes", "no")))

MINPLOT <- min(plot_dat$lci)
MAXPLOT <- max(plot_dat$uci)

plotBetas <- function (v) {
  sub_dat <- plot_dat %>% 
    filter(SNP == v)
  
  res_plot <- ggplot(sub_dat, aes(x = BETA_FLIP, y = strata,
                                  group = term)) +
    geom_pointrange(aes(xmin = lci, xmax = uci,
                        shape = term, color = sex_strata,
                        linetype = sig_lty, alpha = sig_lty),
                    position = position_dodge(width = 0.7),
                    size = 0.4) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_color_manual(values = custom_three_diverge, guide = "none") +
    scale_alpha_manual(values = c(no = 0.7, yes = 1)) +
    scale_x_continuous(limits = c(MINPLOT, MAXPLOT),
                       guide = guide_axis(check.overlap = TRUE)) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  return (res_plot) 
}

tiff(paste0("C:/Users/samvida/Documents/Lindgren Group/Adiposity_Primary_Care/Reports/Manuscript/figures/change_snp_effects/rs429358.tiff"),
     height = 5, width = 5, units = "cm",
     res = 300)
print(plotBetas("rs429358"))
dev.off()

