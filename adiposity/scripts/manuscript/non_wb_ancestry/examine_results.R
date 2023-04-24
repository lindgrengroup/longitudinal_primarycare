# Author: Samvida S. Venkatesh
# Date: 29/09/22

library(tidyverse)
theme_set(theme_bw())

outfile_path <- "" # REDACTED

# colour palette: rose, teal, grey
custom_three_diverge <- c("#D35C79","#009593", "#666666")
names(custom_three_diverge) <- c("F", "M", "sex_comb")

# colour palette: six colours for ancestries
custom_six_diverge <- c("#D35C79", "#D9AB90", "#C7B241", 
                        "#9FBCA4", "#009593", "#7D7D7D")
names(custom_six_diverge) <- c("asian", "black", "chinese",
                               "mixed", "white", "other")

ANC_LEVELS <- c("other", "white", "mixed", "chinese", "black", "asian")

# Read data ----

dat <- read.table("all_ancestries_all_phenos.txt", sep = "\t",
                  header = T, stringsAsFactors = F)

STRATA_LEVELS <- paste0(rep(ANC_LEVELS, each = 3),
                        "_", rep(c("sex_comb", "M", "F"), times = 6))
PTHRESH <- 1E-4

# Self-reported weight change ----

selfrep <- dat %>%
  filter(pheno_tested == "Weight_change_1yr" & visit_compared == 1) 

plot_dat <- selfrep %>%
  mutate(sex_strata = factor(sex_strata, levels = c("F", "M", "sex_comb")),
         ancestry = factor(ancestry, levels = ANC_LEVELS),
         strata = paste0(ancestry, "_", sex_strata),
         sig_lty = factor(ifelse(pval < PTHRESH, "yes", "no"),
                          levels = c("yes", "no")))

plotORs <- function (v) {
  print(v)
  sub_dat <- plot_dat %>% 
    filter(SNP == v & sex_strata == "sex_comb" & !is.na(sig_lty)) %>% 
    select(all_of(c("or", "lci", "uci",
                    "sex_strata", "ancestry", "strata", "sig_lty")))
  missing_anc <- ANC_LEVELS[!ANC_LEVELS %in% unique(sub_dat$ancestry)]

  if (length(missing_anc) > 0) {
    add_dat <- data.frame(or = NA, lci = NA, uci = NA,
                          sex_strata = NA, ancestry = missing_anc,
                          strata = NA, sig_lty = "no")
    sub_dat <- bind_rows(sub_dat, add_dat)
  }
   sub_dat <- sub_dat %>%
    mutate(ancestry = factor(ancestry, levels = ANC_LEVELS))
  
  res_plot <- ggplot(sub_dat, aes(x = or, y = ancestry)) +
    geom_pointrange(aes(xmin = lci, xmax = uci,
                        color = ancestry,
                        linetype = sig_lty, alpha = sig_lty),
                    position = position_dodge(width = 0.7),
                    size = 0.3) +
    geom_vline(xintercept = 1, linetype = 2) +
    scale_color_manual(values = custom_six_diverge, guide = "none") +
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

VARIDS <- unique(plot_dat$SNP)
# Skip rs61955499 because standard erorrs are very small
# VARIDS <- VARIDS[VARIDS != "rs61955499"]
lapply(VARIDS, function (v) {
  tiff(paste0(outfile_path, "/figures/non_wb_ancestry/",
              v, "_selfrep.tiff"),
       height = 3, width = 4, units = "cm",
       res = 300)
  print(plotORs(v))
  dev.off()
})
v <- VARIDS[2]
tiff(paste0(outfile_path, "/figures/non_wb_ancestry/",
            "chr6_26076446", "_selfrep.tiff"),
     height = 3, width = 4, units = "cm",
     res = 300)
print(plotORs(v))
dev.off()

# Linear slope ----

b1_dat <- dat %>%
  filter(term == "b1") 

plot_dat <- b1_dat %>%
  mutate(lci = beta - 1.96*se,
         uci = beta + 1.96*se,
         sex_strata = factor(sex_strata, levels = c("F", "M", "sex_comb")),
         ancestry = factor(ancestry, levels = ANC_LEVELS),
         pheno_tested = factor(pheno_tested, levels = c("Weight", "BMI")),
         strata = factor(paste0(ancestry, "_", pheno_tested),
                         levels = paste0(rep(ANC_LEVELS, each = 2), "_",
                                         rep(c("Weight", "BMI"), times = 6))),
         sig_lty = factor(ifelse(pval < PTHRESH, "yes", "no"),
                          levels = c("yes", "no")))

plotBetas <- function (v) {
  print(v)
  sub_dat <- plot_dat %>% 
    filter(SNP == v & sex_strata == "sex_comb" & !is.na(sig_lty)) %>% 
    select(all_of(c("beta", "lci", "uci",
                    "sex_strata", "ancestry", "pheno_tested", 
                    "strata", "sig_lty")))
  
  missing_anc <- ANC_LEVELS[!ANC_LEVELS %in% unique(sub_dat$ancestry)]
  
  if (length(missing_anc) > 0) {
    add_dat <- data.frame(or = NA, lci = NA, uci = NA,
                          sex_strata = NA, 
                          ancestry = rep(missing_anc, each = 2),
                          pheno_tested = rep(c("BMI", "Weight"), times = length(missing_anc)),
                          strata = NA, sig_lty = "no")
    add_dat$strata <- paste0(add_dat$ancestry, "_", add_dat$pheno_tested)
    sub_dat <- bind_rows(sub_dat, add_dat)
  }
  sub_dat <- sub_dat %>%
    mutate(strata = factor(strata, levels = paste0(rep(ANC_LEVELS, each = 2), "_",
                                                   rep(c("Weight", "BMI"), times = 6))))
  
  res_plot <- ggplot(sub_dat, aes(x = beta, y = strata, group = ancestry)) +
    geom_pointrange(aes(xmin = lci, xmax = uci,
                        shape = pheno_tested, color = ancestry,
                        linetype = sig_lty, alpha = sig_lty),
                    position = position_dodge(width = 0.7),
                    size = 0.3) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_color_manual(values = custom_six_diverge, guide = "none") +
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

VARIDS <- unique(plot_dat$SNP)
# Skip rs61955499 because sample sizes are very small
# VARIDS <- VARIDS[VARIDS != "rs61955499"]
lapply(VARIDS, function (v) {
  tiff(paste0(outfile_path, "/figures/non_wb_ancestry/",
              v, "_b1.tiff"),
       height = 3, width = 4, units = "cm",
       res = 300)
  print(plotBetas(v))
  dev.off()
})
v <- VARIDS[2]
tiff(paste0(outfile_path, "/figures/non_wb_ancestry/",
            "chr6_26076446", "_b1.tiff"),
     height = 3, width = 4, units = "cm",
     res = 300)
print(plotBetas(v))
dev.off()

# p(k1) ----

k1_dat <- dat %>%
  filter(term == "k1") 

plot_dat <- k1_dat %>%
  mutate(sex_strata = factor(sex_strata, levels = c("F", "M", "sex_comb")),
         ancestry = factor(ancestry, levels = ANC_LEVELS),
         pheno_tested = factor(pheno_tested, levels = c("Weight", "BMI")),
         strata = factor(paste0(ancestry, "_", pheno_tested),
                         levels = paste0(rep(ANC_LEVELS, each = 2), "_",
                                         rep(c("Weight", "BMI"), times = 6))),
         sig_lty = factor(ifelse(pval < PTHRESH, "yes", "no"),
                          levels = c("yes", "no")))

plotORs <- function (v) {
  print(v)
  sub_dat <- plot_dat %>% 
    filter(SNP == v & sex_strata == "sex_comb" & !is.na(sig_lty)) %>% 
    select(all_of(c("or", "lci", "uci",
                    "sex_strata", "ancestry", "pheno_tested", 
                    "strata", "sig_lty")))
  
  missing_anc <- ANC_LEVELS[!ANC_LEVELS %in% unique(sub_dat$ancestry)]
  
  if (length(missing_anc) > 0) {
    add_dat <- data.frame(or = NA, lci = NA, uci = NA,
                          sex_strata = NA, 
                          ancestry = rep(missing_anc, each = 2),
                          pheno_tested = rep(c("BMI", "Weight"), times = length(missing_anc)),
                          strata = NA, sig_lty = "no")
    add_dat$strata <- paste0(add_dat$ancestry, "_", add_dat$pheno_tested)
    sub_dat <- bind_rows(sub_dat, add_dat)
  }
  sub_dat <- sub_dat %>%
    mutate(strata = factor(strata, levels = paste0(rep(ANC_LEVELS, each = 2), "_",
                                                   rep(c("Weight", "BMI"), times = 6))))
  
  
  res_plot <- ggplot(sub_dat, aes(x = or, y = strata, group = ancestry)) +
    geom_pointrange(aes(xmin = lci, xmax = uci,
                        shape = pheno_tested, color = ancestry,
                        linetype = sig_lty, alpha = sig_lty),
                    position = position_dodge(width = 0.7),
                    size = 0.3) +
    geom_vline(xintercept = 1, linetype = 2) +
    scale_color_manual(values = custom_six_diverge, guide = "none") +
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

VARIDS <- unique(plot_dat$SNP)
# Skip rs61955499 because sample sizes are very small
# VARIDS <- VARIDS[VARIDS != "rs61955499"]
lapply(VARIDS, function (v) {
  tiff(paste0(outfile_path, "/figures/non_wb_ancestry/",
              v, "_k1.tiff"),
       height = 3, width = 4, units = "cm",
       res = 300)
  print(plotORs(v))
  dev.off()
})
v <- VARIDS[2]
tiff(paste0(outfile_path, "/figures/non_wb_ancestry/",
            "chr6_26076446", "_k1.tiff"),
     height = 3, width = 4, units = "cm",
     res = 300)
print(plotORs(v))
dev.off()

# Write for tables ----

to_write <- dat %>%
  filter(visit_compared == 1 | is.na(visit_compared)) %>%
  filter(is.na(term) | term %in% c("b1", "k1")) %>%
  mutate(beta_se = paste0(signif(beta, 3), " (", signif(se, 3), ")"),
         or_ci = paste0(signif(or, 3), " (", signif(lci, 3), " - ", signif(uci, 3), ")"),
         write_p = signif(pval, 3),
         strata = factor(paste0(pheno_tested, "_", term),
                         levels = c(paste0(rep(c("BMI", "Weight"), each = 2), "_",
                                           rep(c("b1", "k1"), times = 2)),
                                    "Weight_change_1yr_NA")),
         SNP = factor(SNP, levels = c("rs9467663", "6:26076446_GA_G",
                                      "rs11778922", "rs61955499", "rs12953815",
                                      "rs429358"))) %>%
  mutate(beta_se = ifelse(strata %in% c("BMI_b1", "Weight_b1"), beta_se, NA)) %>%
  arrange(SNP, ancestry, strata) %>%
  select(all_of(c("SNP", "ancestry", "strata", "sex_strata",
                  "beta_se", "or_ci", "write_p", "sample_size")))

to_write <- to_write %>% 
  pivot_wider(id_cols = c(SNP, ancestry, strata),
              names_from = sex_strata,
              values_from = c(beta_se, or_ci, write_p, sample_size))

to_write <- to_write[, c("SNP", "ancestry", "strata", 
                         "beta_se_sex_comb", "or_ci_sex_comb", "write_p_sex_comb", "sample_size_sex_comb",
                         "beta_se_F", "or_ci_F", "write_p_F", "sample_size_F",
                         "beta_se_M", "or_ci_M", "write_p_M", "sample_size_M")]

write.table(to_write, "all_results_all_phenos_for_tables.txt",
            sep = "\t", row.names = F, quote = F)
