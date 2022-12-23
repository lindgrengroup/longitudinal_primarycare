# Author: Samvida S. Venkatesh
# Date: 29/09/22

library(tidyverse)
theme_set(theme_bw())

# colour palette: rose, teal, grey
custom_three_diverge <- c("#D35C79","#009593", "#666666")
names(custom_three_diverge) <- c("F", "M", "sex_comb")

# Read data ----

dat <- read.table("rs429358_replication_by_dementia_status.txt", sep = "\t",
                  header = T, stringsAsFactors = F)
PTHRESH <- 0.01

# Self-reported weight change ----

selfrep <- dat %>%
  filter(pheno_tested == "Weight_change_1yr" & visit_compares == 1) 

plot_dat <- selfrep %>%
  mutate(sex_strata = factor(sex_strata, levels = c("F", "M", "sex_comb")),
         strata = factor(paste0(sex_strata, "_", status),
                         levels = c("sex_comb_no dementia", "sex_comb_all",
                                    "M_no dementia", "M_all",
                                    "F_no dementia", "F_all")),
         sig_lty = factor(ifelse(pval < PTHRESH, "yes", "no"),
                          levels = c("yes", "no")))

selfrep_plot <- ggplot(plot_dat, aes(x = or, y = strata)) +
  geom_pointrange(aes(xmin = lci, xmax = uci,
                      color = sex_strata, shape = status,
                      linetype = sig_lty, alpha = sig_lty),
                  position = position_dodge(width = 0.7),
                  size = 0.3) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_color_manual(values = custom_three_diverge, guide = "none") +
  scale_alpha_manual(values = c(no = 0.7, yes = 1)) +
  scale_linetype_manual(values = c(no = 2, yes = 1)) +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

tiff("C:/Users/samvida/Documents/Lindgren Group/Adiposity_Primary_Care/Reports/Manuscript/figures/rs429358_dementia/selfrep_wtchg.tiff",
     height = 5, width = 5, units = "cm",
     res = 300)
print(selfrep_plot)
dev.off()

# Linear slope ----

b1_dat <- dat %>%
  filter(term == "b1") 

plot_dat <- b1_dat %>%
  mutate(lci = beta - 1.96*se,
         uci = beta + 1.96*se,
         sex_strata = factor(sex_strata, levels = c("F", "M", "sex_comb")),
         strata = factor(paste0(sex_strata, "_", status),
                         levels = c("sex_comb_no dementia", "sex_comb_all",
                                    "M_no dementia", "M_all",
                                    "F_no dementia", "F_all")),
         sig_lty = factor(ifelse(pval < PTHRESH, "yes", "no"),
                          levels = c("yes", "no")))

plotBetas <- function (p) {
  sub_dat <- plot_dat %>% 
    filter(pheno_tested == p) 
  
  res_plot <- ggplot(sub_dat, aes(x = beta, y = strata, 
                                  group = sex_strata)) +
    geom_pointrange(aes(xmin = lci, xmax = uci,
                        shape = status, color = sex_strata,
                        linetype = sig_lty, alpha = sig_lty),
                    position = position_dodge(width = 0.7),
                    size = 0.3) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_color_manual(values = custom_three_diverge, guide = "none") +
    scale_alpha_manual(values = c(no = 0.7, yes = 1)) +
    scale_linetype_manual(values = c(no = 2, yes = 1)) +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  return (res_plot) 
}

PHENOTYPES <- c("BMI", "Weight")
lapply(PHENOTYPES, function (p) {
  tiff(paste0("C:/Users/samvida/Documents/Lindgren Group/Adiposity_Primary_Care/Reports/Manuscript/figures/rs429358_dementia/",
              p, "_b1.tiff"),
       height = 5, width = 5, units = "cm",
       res = 300)
  print(plotBetas(p))
  dev.off()
})

# p(k1) ----

k1_dat <- dat %>%
  filter(term == "k1") 

plot_dat <- k1_dat %>%
  mutate(sex_strata = factor(sex_strata, levels = c("F", "M", "sex_comb")),
         strata = factor(paste0(sex_strata, "_", status),
                         levels = c("sex_comb_no dementia", "sex_comb_all",
                                    "M_no dementia", "M_all",
                                    "F_no dementia", "F_all")),
         sig_lty = factor(ifelse(pval < PTHRESH, "yes", "no"),
                          levels = c("yes", "no")))


plotORs <- function (p) {
  sub_dat <- plot_dat %>% 
    filter(pheno_tested == p) 
  
  res_plot <- ggplot(sub_dat, aes(x = or, y = strata, 
                                  group = sex_strata)) +
    geom_pointrange(aes(xmin = lci, xmax = uci,
                        shape = status, color = sex_strata,
                        linetype = sig_lty, alpha = sig_lty),
                    position = position_dodge(width = 0.7),
                    size = 0.3) +
    geom_vline(xintercept = 1, linetype = 2) +
    scale_color_manual(values = custom_three_diverge, guide = "none") +
    scale_alpha_manual(values = c(no = 0.7, yes = 1)) +
    scale_linetype_manual(values = c(no = 2, yes = 1)) +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  return (res_plot) 
}

PHENOTYPES <- c("BMI", "Weight")
lapply(PHENOTYPES, function (p) {
  tiff(paste0("C:/Users/samvida/Documents/Lindgren Group/Adiposity_Primary_Care/Reports/Manuscript/figures/rs429358_dementia/",
              p, "_k1.tiff"),
       height = 5, width = 5, units = "cm",
       res = 300)
  print(plotORs(p))
  dev.off()
})

