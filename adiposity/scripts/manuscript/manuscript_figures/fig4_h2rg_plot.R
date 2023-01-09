# Author: Samvida S. Venkatesh
# Date: 08/09/2022

library(tidyverse)
theme_set(theme_bw())

# colour palette: rose, teal, grey
custom_three_diverge <- c("#D35C79","#009593", "#666666")
names(custom_three_diverge) <- c("F", "M", "sex_comb")

# Read data ----

plot_dir <- "" # REDACTED

dat <- read.table("heritability_rg.txt", sep = "\t", header = T,
                  stringsAsFactors = F)

dat <- dat %>% filter(term %in% c("b0", "b1", "k1"))

# Heritability forest plot ----

h2_dat <- dat %>%
  mutate(sex_strata = factor(sex_strata, levels = c("F", "M", "sex_comb")),
         strata = factor(paste0(phenotype, "_", sex_strata), levels = c("Weight_sex_comb", "Weight_M", "Weight_F",
                                                                        "BMI_sex_comb", "BMI_M", "BMI_F")),
         uci_h2 = pmin(beta_h2 + 1.96*se_h2, 1),
         lci_h2 = pmax(beta_h2 - 1.96*se_h2, 0))

h2_plot <- ggplot(h2_dat, aes(x = beta_h2, y = strata)) +
  facet_wrap(~term, nrow = 3, ncol = 1, scales = "free") +
  geom_pointrange(aes(xmin = lci_h2, xmax = uci_h2,
                      shape = phenotype, color = sex_strata),
                  size = 0.3,
                  position = position_dodge(width = 0.5)) +
  scale_color_manual(values = custom_three_diverge, guide = "none") +
  scale_x_continuous(limits = c(0, 0.35),
                     guide = guide_axis(check.overlap = TRUE)) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        strip.background = element_rect(size = 0.7),
        strip.text.x = element_text(size = 8))

tiff(paste0(plot_dir, "/h2rg/heritability.tiff"),
     height = 10, width = 5, units = "cm",
     res = 300)
print(h2_plot)
dev.off()

h2_sexhet <- h2_dat %>%
  pivot_wider(id_cols = c(phenotype, term), 
              names_from = sex_strata, 
              values_from = c(beta_h2, se_h2)) %>%
  mutate(zstat_sexhet = abs((beta_h2_F - beta_h2_M) / sqrt(se_h2_F^2 + se_h2_M^2)),
         pval_sexhet = 1 - pnorm(zstat_sexhet))


rg_dat <- dat %>%
  mutate(sex_strata = factor(sex_strata, levels = c("F", "M", "sex_comb")),
         strata = factor(paste0(phenotype, "_", sex_strata), levels = c("Weight_sex_comb", "Weight_M", "Weight_F",
                                                                        "BMI_sex_comb", "BMI_M", "BMI_F")),
         uci_rg = pmin(beta_rg + 1.96*se_rg, 1),
         lci_rg = pmax(beta_rg - 1.96*se_rg, 0))

rg_plot <- ggplot(rg_dat, aes(x = beta_rg, y = strata)) +
  facet_wrap(~term, nrow = 3, ncol = 1, scales = "free") +
  geom_pointrange(aes(xmin = lci_rg, xmax = uci_rg,
                      shape = phenotype, color = sex_strata),
                  size = 0.3,
                  position = position_dodge(width = 0.5)) +
  scale_color_manual(values = custom_three_diverge, guide = "none") +
  scale_x_continuous(limits = c(0,1),
                     guide = guide_axis(check.overlap = TRUE)) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        strip.background = element_rect(size = 0.7),
        strip.text.x = element_text(size = 8))


tiff(paste0(plot_dir, "/h2rg/rg.tiff"),
     height = 10, width = 5, units = "cm",
     res = 300)
print(rg_plot)
dev.off()

rg_sexhet <- rg_dat %>%
  pivot_wider(id_cols = c(phenotype, term), 
              names_from = sex_strata, 
              values_from = c(beta_rg, se_rg)) %>%
  mutate(zstat_sexhet = abs((beta_rg_F - beta_rg_M) / sqrt(se_rg_F^2 + se_rg_M^2)),
         pval_sexhet = 1 - pnorm(zstat_sexhet))

