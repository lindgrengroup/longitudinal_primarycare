# Author: Samvida S. Venkatesh
# Date: 31/08/22

library(tidyverse)
library(lme4)
library(broom)
library(foreign)
library(MASS)
library(ggpubr)
theme_set(theme_bw())

# Read in arguments ----

plotdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/plots/"
resdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/replication_results/"

# For 3 colours that are sequential teal
custom_teal_sequential <- c("#005958", "#009593", "#66BFBE")
# For 3 colours that are sequential rose
custom_rose_sequential <- c("#7E3748", "#D35C79", "#E49DAE")
# For 3 colours that are sequential black/grey
custom_black_sequential <- c("#000000", "#666666", "#CCCCCC")
# colour palette: rose, yellow, teal
custom_three_diverge <- c("#D35C79","#C7B241", "#009593")
names(custom_three_diverge) <- c("Gain", "No change", "Loss")

# Read data ----

main_longit_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/data/main_data_adipo_change.rds")
QUANT_PHENOS <- c("BMI", "Weight", "WC", "WHR",
                  "WCadjBMI", "WHRadjBMI")
CAT_PHENOS <- "Weight_change_1yr"

gp_ids <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/data/all_ids_in_discovery_gp_ukb_dat.txt",
                     sep = "\t", header = F, stringsAsFactors = F)$V1
gp_ids <- as.character(gp_ids)

vars_to_replicate <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/data/lead_snps_to_replicate.txt",
                                sep = "\t", header = T, stringsAsFactors = F)
VARIDS <- vars_to_replicate$SNP
# For variants without rsids, build variant name
rename_vars <- grep("^chr", VARIDS)
VARIDS[rename_vars] <- paste0(vars_to_replicate$SNP[rename_vars], "_",
                              vars_to_replicate$Tested_Allele[rename_vars], "_", 
                              vars_to_replicate$Other_Allele[rename_vars])
VARIDS <- gsub("^chr", "", VARIDS)

# Genotypes / dosages at variants of interest
var_dosages <- lapply(VARIDS, function (varid) {
  res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sample_variant_counts/",
                           varid, "_dosages.txt"),
                    sep = " ", header = T, stringsAsFactors = F)
  # Remove first row, which contains info on type of column and columns 
  # 2, 3, 4 (ID repeat, missingness, sex)
  res <- res[-1, c(1, 5)]
  colnames(res) <- c("eid", varid)
  return (res)
})
names(var_dosages) <- VARIDS

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

MOD_COVARS <- c("sex", "baseline_age", "age_sq", "FUyrs",
                "year_of_birth", "data_provider", 
                paste0("PC", 1:21)) # add BMI as covariate for categorical tests

general_covars <- general_covars %>% 
  dplyr::select(any_of(c("eid", MOD_COVARS))) %>%
  mutate(sex = factor(sex), year_of_birth = as.numeric(year_of_birth))

SEX_STRATA <- c("F", "M", "sex_comb")

# Wrangle data ----

# Convert dosages to 0/1/2 genotypes based on threshold
DOSAGE_THRESHOLD_0 <- 0.5
DOSAGE_THRESHOLD_2 <- 1.5
var_dosages_hardcall <- lapply(VARIDS, function (varid) {
  res <- var_dosages[[varid]] %>% 
    mutate(genotype = ifelse(!!as.symbol(varid) < DOSAGE_THRESHOLD_0, "0",
                             ifelse(!!as.symbol(varid) > DOSAGE_THRESHOLD_2, 
                                    "2", "1")),
           eid = as.character(eid)) %>%
    dplyr::select(all_of(c("eid", "genotype")))
  return (res)
})
names(var_dosages_hardcall) <- VARIDS

# Remove IDs with GP data from main longitudinal data for replication
# Tally visit number and changes from first visit needed for model
# Add in model covariates

dat_for_replication <- lapply(c(QUANT_PHENOS, CAT_PHENOS), function (p) {
  dat <- main_longit_dat[[p]] %>% 
    filter(!eid %in% gp_ids) %>%
    group_by(eid) %>% 
    mutate(visit = seq(n()),
           baseline_age = first(age_event),
           age_sq = baseline_age^2,
           age_event_sq = age_event^2,
           FUyrs = age_event - first(age_event)) %>%
    left_join(general_covars)
  
  if (p %in% QUANT_PHENOS) 
    dat <- dat %>% mutate(value_diff = value - first(value))
  else if (p %in% CAT_PHENOS)
    dat <- dat %>% mutate(value = factor(value, 
                                         levels = c("Loss", "No change", "Gain"),
                                         ordered = T))
  
  dat <- dat[complete.cases(dat), ]
  return (dat)
})
names(dat_for_replication) <- c(QUANT_PHENOS, CAT_PHENOS)

# Function to add variant dosage
addVarDosage <- function (dat, varid) {
  res <- dat
  var_dat <- var_dosages[[varid]]
  res$dosage <- var_dat[match(res$eid, var_dat$eid), varid]
  res$dosage <- as.numeric(res$dosage)
  res <- res[!is.na(res$dosage), ]
  return (res)
}

# Function to add hard-called genotype
addGenoGroup <- function (dat, varid) {
  res <- dat
  var_dat <- var_dosages_hardcall[[varid]]
  res$genotype <- var_dat$genotype[match(res$eid, var_dat$eid)]
  res$genotype <- factor(res$genotype, levels = c("0", "1", "2"))
  res <- res[!is.na(res$genotype), ]
  return (res)
}

# Testing functions ----

getLMMTerm <- function (dat, mod_cvs) {
  to_model <- dat %>%
    group_by(eid) %>%
    arrange(age_event, .by_group = T) %>%
    mutate(t = age_event - first(age_event))
  
  mod_form <- paste0("value ~ t + ",
                     paste0(mod_cvs, collapse = " + "),
                     " + (t | eid)")
  # Run model
  lmod <- lmer(formula(mod_form), data = to_model, REML = F)
  # Get BLUPs 
  blups <- coef(lmod)$eid[, c("(Intercept)", "t")]
  colnames(blups) <- c("intercept", "slope")
  
  # Adjust slope term for intercept and RINT
  modelled_slope <- lm(slope ~ intercept, blups)
  
  # RINT
  trait_to_rint <- resid(modelled_slope)
  rinted_trait <- qnorm((rank(trait_to_rint) - 0.5) / sum(!is.na(trait_to_rint)))
  
  res <- data.frame(eid = rownames(blups),
                    lmm_slope_adj_int = rinted_trait)
  res <- inner_join(res, dat, by = "eid")
  return (res)
}

quantTest <- function (dat, ss) {
  covars_include <- MOD_COVARS
  if (ss != "sex_comb") {
    dat <- dat %>% filter(sex == ss)
    covars_include <- covars_include[-which(covars_include == "sex")]
  }
  if (length(unique(dat$data_provider)) < 2)
    covars_include <- covars_include[-which(covars_include == "data_provider")]
  
  # Get term to use for association testing
  # i.e. RINTed adjusted LMM slope
  dat_final <- getLMMTerm(dat, covars_include)
  
  modeled_dat <- lm(lmm_slope_adj_int ~ dosage, data = dat_final)
  
  print_res <- tidy(modeled_dat) %>%
    filter(term == "dosage") %>%
    rename(beta = estimate, se = std.error, tstat = statistic, pval = p.value) %>%
    dplyr::select(all_of(c("beta", "se", "tstat", "pval"))) %>%
    mutate(sample_size = nrow(dat))
  
  return (print_res)
}

catTest <- function (dat, ss) {
  # Get the correct covariates for adjustment
  covars_include <- c("sex", "age_event", "age_event_sq",
                      "year_of_birth", "data_provider", "BMI",
                      paste0("PC", 1:21))
  if (ss != "sex_comb") {
    dat <- dat %>% filter(sex == ss)
    covars_include <- covars_include[-which(covars_include == "sex")]
  }
  if (length(unique(dat$data_provider)) < 2)
    covars_include <- covars_include[-which(covars_include == "data_provider")]
  
  mod_formula <- paste0("value ~ dosage + ",
                        paste0(covars_include, collapse = " + "))
  
  modeled_dat <- polr(formula(mod_formula), data = dat, Hess = T)
  print_res <- tidy(modeled_dat) %>%
    filter(term == "dosage") %>%
    rename(beta = estimate, se = std.error, tstat = statistic) %>%
    dplyr::select(all_of(c("beta", "se", "tstat"))) %>%
    mutate(sample_size = nrow(dat))
  return (print_res)
}

# Plotting functions ----

quantPlot <- function (dat, qp, varid) {
  
  per_sex <- lapply(SEX_STRATA, function (ss) {
    # Adjust value-diff for covariates
    covars_adj <- MOD_COVARS
    if (ss != "sex_comb") {
      dat <- dat %>% filter(sex == ss)
      covars_adj <- covars_adj[-which(covars_adj == "sex")]
    }
    if (length(unique(dat$data_provider)) < 2)
      covars_adj <- covars_adj[-which(covars_adj == "data_provider")]
    
    mod_formula <- paste0("value_diff ~ ",
                          paste0(covars_adj, collapse = " + "))
    modeled_dat <- lm(formula(mod_formula), data = dat)
    dat$value_diff_adj <- resid(modeled_dat)
    
    summ_dat <- dat %>%
      filter(visit != 1) %>%
      group_by(visit, genotype) %>%
      summarise(mean_se = mean_se(value_diff_adj, 1.96)) %>%
      unnest(mean_se) %>%
      mutate(visit = factor(visit, levels = c("2", "3")),
             dosage_gp = factor(genotype, levels = c("0", "1", "2")))
    
    if (ss == "F") use_col_palette <- custom_rose_sequential
    else if (ss == "M") use_col_palette <- custom_teal_sequential
    else use_col_palette <- custom_black_sequential
    
    resplot <- ggplot(summ_dat, aes(x = y, y = visit, 
                                    fill = genotype, color = genotype)) +
      geom_pointrange(aes(xmin = ymin, xmax = ymax),
                      position = position_dodge(width = 0.7)) + 
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_color_manual(values = use_col_palette) +
      scale_fill_manual(values = use_col_palette) +
      scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
      labs(x = paste0("Change in ", qp, " from visit 1"), 
           title = paste0(varid, ": ", ss))
    return (resplot)
  })
  
  final_plot <- ggarrange(plotlist = per_sex,
                          nrow = 3, ncol = 1, common.legend = T)
  return (final_plot)
}

catPlot <- function (dat) {
  summ_dat <- dat %>%
    group_by(visit, genotype) %>%
    count(value) %>%
    ungroup() %>% group_by(visit, genotype) %>%
    mutate(prop = n/sum(n)) %>% 
    filter(value %in% c("Gain", "No change", "Loss")) %>%
    mutate(value = factor(value, 
                          levels = c("Gain", "No change", "Loss")),
           visit = factor(visit, levels = c("1", "2", "3")),
           genotype = factor(genotype, levels = c("0", "1", "2")))
  
  resplot <- ggplot(summ_dat, aes(x = genotype, y = n,
                                  fill = value, color = value)) +
    facet_wrap(~visit, ncol = 3) +
    geom_bar(position = "fill", stat = "identity") + 
    scale_color_manual(values = custom_three_diverge) +
    scale_fill_manual(values = custom_three_diverge) +
    scale_y_continuous(guide = guide_axis(check.overlap = TRUE))
  
  return (resplot)
}

# Apply testing and plotting ----

adipo_change_replication <- lapply(VARIDS, function (varid) {
  cat(paste0("Running SNP: ", varid, "\n"))
  dir.create(paste0(plotdir, varid))
  # Loop through the various adiposity traits
  
  quant_res_tables <- lapply(QUANT_PHENOS, function (qp) {
    # Create full data with dosage and genotype info
    model_dat <- dat_for_replication[[qp]]
    model_dat <- addVarDosage(model_dat, varid)
    model_dat <- addGenoGroup(model_dat, varid)
    
    # Run regressions per sex
    regression_res <- lapply(SEX_STRATA, function (ss) {
      res <- quantTest(model_dat, ss) %>%
        mutate(sex_strata = ss)
      return (res)
    })
    regression_res <- bind_rows(regression_res) %>%
      mutate(pheno_tested = qp)
    
    # Plot change in adiposity
    png(paste0(plotdir, varid, "/", varid, "_", qp, ".png"))
    print(quantPlot(model_dat, qp, varid))
    dev.off()
    
    return (regression_res)
  })
  quant_res_tables <- bind_rows(quant_res_tables)
  
  # Loop through categorical traits (only one for now)
  cat_res_tables <- lapply(CAT_PHENOS, function (cp) {
    # Create full data with dosage and genotype info
    model_dat <- dat_for_replication[[cp]]
    model_dat <- addVarDosage(model_dat, varid)
    model_dat <- addGenoGroup(model_dat, varid)
    
    # Run regressions and plots
    regression_res <- lapply(SEX_STRATA, function (ss) {
      per_visit <- lapply(c(1:3), function (vc) {
        sub_dat <- model_dat %>% filter(visit == vc)
        res <- catTest(sub_dat, ss) %>%
          mutate(visit_compared = vc)
      })
      per_visit <- bind_rows(per_visit) %>%
        mutate(sex_strata = ss)
      
      # Plot
      plot_dat <- model_dat
      if (ss != "sex_comb") plot_dat <- plot_dat %>% filter(sex == ss)
      png(paste0(plotdir, varid, "/", varid, "_", cp, "_", ss, ".png"))
      print(catPlot(plot_dat) + 
              labs(y = cp, title = paste0(varid, ": ", ss)))
      dev.off()
      
      return (per_visit)
    })
    regression_res <- bind_rows(regression_res) %>%
      mutate(pheno_tested = cp)
    return (regression_res)
  })
  cat_res_tables <- bind_rows(cat_res_tables) %>%
    mutate(OR = exp(beta), 
           lci = exp(beta - (1.96*se)), uci = exp(beta + (1.96*se)),
           pval = pt(tstat, df = sample_size))
  
  full_res_tables <- bind_rows(quant_res_tables, cat_res_tables) %>%
    mutate(SNP = varid)
  return (full_res_tables)
})
adipo_change_replication <- bind_rows(adipo_change_replication)

write.table(adipo_change_replication, 
            paste0(resdir, "adipo_change_replication.txt"),
            sep = "\t", quote = F, row.names = F)

# Explore results ----

dat <- read.table("adipo_change_replication_lmm_slopes.txt", sep = "\t",
                  header = T, stringsAsFactors = F)

# Subset to phenotypes of interest at visits with most power (sample size)

selfrep_wtchg <- dat %>% 
  filter(pheno_tested == "Weight_change_1yr" & visit_compared == 1)

sig_wtchg_snps <- selfrep_wtchg %>%
  filter(pval <= 0.005) %>%
  mutate(perc_change = signif((1 - OR)*100, 3))

# Subset to abdominal obesity
abdo_obesity <- dat %>%
  filter(pheno_tested %in% c("WC", "WCadjBMI", "WHR", "WHRadjBMI")) 

sig_abdo_obesity <- abdo_obesity %>%
  filter(pval <= 0.05/4)

# printing to table
to_write <- bind_rows(selfrep_wtchg, abdo_obesity) %>%
  filter(grepl("^rs", SNP)) %>%
  dplyr::select(all_of(c("beta", "se", "pval", "sex_strata", "pheno_tested",
                         "OR", "lci", "uci", "SNP"))) %>%
  mutate(beta_print = paste0(signif(beta, 3), " (", signif(se, 3), ")"),
         OR_print = paste0(signif(OR, 3), " (", signif(lci, 3), " - ", signif(uci, 3), ")"),
         pval_print = signif(pval, 3)) %>%
  dplyr::select(all_of(c("SNP", "pheno_tested", "sex_strata", "beta_print", "OR_print", "pval_print"))) %>%
  arrange(SNP, sex_strata, factor(pheno_tested, levels = c("WC", "WCadjBMI",
                                                           "WHR", "WHRadjBMI",
                                                           "Weight_change_1yr")))
write.table(to_write, "for_results_tables.txt",
            sep = "\t", quote = F, row.names = F)

# Plot results ----

theme_set(theme_bw())

# colour palette: rose, teal, grey
custom_three_diverge <- c("#D35C79","#009593", "#666666")
names(custom_three_diverge) <- c("F", "M", "sex_comb")

plot_dat <- abdo_obesity %>%
  filter(grepl("^rs", SNP)) %>%
  mutate(sex_strata = factor(sex_strata, levels = c("F", "M", "sex_comb")),
         pheno_tested = factor(pheno_tested, levels = c("WHRadjBMI", "WHR",
                                                        "WCadjBMI", "WC")),
         bmi_adj = factor(ifelse(grepl("adjBMI", pheno_tested), "yes", "no")),
         strata = factor(paste0(pheno_tested, "_", sex_strata), 
                         levels = c("WHRadjBMI_sex_comb", "WHR_sex_comb",
                                    "WHRadjBMI_M", "WHR_M",
                                    "WHRadjBMI_F", "WHR_F",
                                    "WCadjBMI_sex_comb", "WC_sex_comb",
                                    "WCadjBMI_M", "WC_M",
                                    "WCadjBMI_F", "WC_F")),
         uci = beta + 1.96*se,
         lci = beta - 1.96*se,
         sig_lty = factor(ifelse(pval < 0.05/4, "yes", "no"),
                          levels = c("yes", "no")))

MINPLOT <- min(plot_dat$lci)
MAXPLOT <- max(plot_dat$uci)

VARIDS <- unique(plot_dat$SNP)

plotBetas <- function (v) {
  sub_dat <- plot_dat %>% 
    filter(SNP == v)
  
  res_plot <- ggplot(sub_dat, aes(x = beta, y = strata,
                                  group = pheno_tested)) +
    geom_pointrange(aes(xmin = lci, xmax = uci,
                        shape = bmi_adj, color = sex_strata,
                        linetype = sig_lty, alpha = sig_lty),
                    position = position_dodge(width = 0.7),
                    size = 0.4) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_color_manual(values = custom_three_diverge, guide = "none") +
    scale_alpha_manual(values = c(no = 0.7, yes = 1)) +
    scale_linetype_manual(values = c(no = 2, yes = 1)) +
    scale_x_continuous(limits = c(MINPLOT, MAXPLOT),
                       guide = guide_axis(check.overlap = TRUE)) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  return (res_plot) 
}

lapply(VARIDS, function (v) {
  tiff(paste0("C:/Users/samvida/Documents/Lindgren Group/Adiposity_Primary_Care/Reports/Manuscript/figures/abdominal_obesity_effects/",
              v, ".tiff"),
       height = 4, width = 5, units = "cm",
       res = 300)
  print(plotBetas(v))
  dev.off()
})

# Plot concordant effects on BMI and WHRadjBMI (scatter with lines for C.I.s) ----

for_plot <- dat %>%
  filter(pheno_tested %in% c("WCadjBMI", "WHRadjBMI", "BMI") & sex_strata == "sex_comb") %>%
  select(all_of(c("beta", "se", "pval", "pheno_tested", "SNP"))) %>%
  mutate(lci = beta - 1.96*se,
         uci = beta + 1.96*se,
         sig_lty = factor(ifelse(pval < 2E-3, "yes", "no"),
                                    levels = c("yes", "no")))

for_plot <- for_plot %>%
  pivot_wider(id_cols = SNP,
              names_from = pheno_tested,
              values_from = c(beta, se, pval, lci, uci, sig_lty))

# Align to the BMI-increasing allele
for_plot <- for_plot %>%
  mutate(flip = beta_BMI < 0,
         og_beta_BMI = beta_BMI,
         beta_BMI = ifelse(flip, -og_beta_BMI, og_beta_BMI),
         lci_BMI = ifelse(flip, -og_beta_BMI - 1.96*se_BMI, lci_BMI),
         uci_BMI = ifelse(flip, -og_beta_BMI + 1.96*se_BMI, uci_BMI),
         og_beta_WHRadjBMI = beta_WHRadjBMI,
         beta_WHRadjBMI = ifelse(flip, -og_beta_WHRadjBMI, og_beta_WHRadjBMI),
         lci_WHRadjBMI = ifelse(flip, -og_beta_WHRadjBMI - 1.96*se_WHRadjBMI, lci_WHRadjBMI),
         uci_WHRadjBMI = ifelse(flip, -og_beta_WHRadjBMI + 1.96*se_WHRadjBMI, uci_WHRadjBMI),
         og_beta_WCadjBMI = beta_WCadjBMI,
         beta_WCadjBMI = ifelse(flip, -og_beta_WCadjBMI, og_beta_WCadjBMI),
         lci_WCadjBMI = ifelse(flip, -og_beta_WCadjBMI - 1.96*se_WCadjBMI, lci_WCadjBMI),
         uci_WCadjBMI = ifelse(flip, -og_beta_WCadjBMI + 1.96*se_WCadjBMI, uci_WCadjBMI))

minplot <- min(c(for_plot$lci_BMI, for_plot$lci_WHRadjBMI))
maxplot <- max(c(for_plot$uci_BMI, for_plot$uci_WHRadjBMI))

whradjbmi_plot <- ggplot(for_plot, aes(x = beta_BMI, y = beta_WHRadjBMI)) +
  geom_vline(xintercept = 0, linetype = 1, color = "grey",
              size = 0.3) +
  geom_hline(yintercept = 0, linetype = 1, color = "grey",
             size = 0.3) +
  geom_pointrange(aes(xmin = lci_BMI, xmax = uci_BMI,
                      linetype = sig_lty_WHRadjBMI, alpha = sig_lty_WHRadjBMI),
                  size = 0.3, fatten = 0.5) +
  geom_pointrange(aes(ymin = lci_WHRadjBMI, ymax = uci_WHRadjBMI,
                      linetype = sig_lty_WHRadjBMI, alpha = sig_lty_WHRadjBMI),
                  size = 0.3, fatten = 0.5) +
  scale_alpha_manual(values = c(no = 0.5, yes = 1)) +
  scale_linetype_manual(values = c(no = 2, yes = 1)) +
  scale_x_continuous(limits = c(minplot, maxplot),
                     guide = guide_axis(check.overlap = TRUE)) +
  scale_y_continuous(limits = c(minplot, maxplot),
                     guide = guide_axis(check.overlap = TRUE)) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

tiff("C:/Users/samvida/Documents/Lindgren Group/Adiposity_Primary_Care/Reports/Manuscript/figures/abdominal_obesity_effects/WHRadjBMI_vs_BMI.tiff",
     height = 5, width = 5, units = "cm",
     res = 300)
print(whradjbmi_plot)
dev.off()

minplot <- min(c(for_plot$lci_BMI, for_plot$lci_WCadjBMI))
maxplot <- max(c(for_plot$uci_BMI, for_plot$uci_WCadjBMI))

wcadjbmi_plot <- ggplot(for_plot, aes(x = beta_BMI, y = beta_WCadjBMI)) +
  geom_vline(xintercept = 0, linetype = 1, color = "grey",
             size = 0.3) +
  geom_hline(yintercept = 0, linetype = 1, color = "grey",
             size = 0.3) +
  geom_pointrange(aes(xmin = lci_BMI, xmax = uci_BMI,
                      linetype = sig_lty_WCadjBMI, alpha = sig_lty_WCadjBMI),
                  size = 0.3, fatten = 0.5) +
  geom_pointrange(aes(ymin = lci_WCadjBMI, ymax = uci_WCadjBMI,
                      linetype = sig_lty_WCadjBMI, alpha = sig_lty_WCadjBMI),
                  size = 0.3, fatten = 0.5) +
  scale_alpha_manual(values = c(no = 0.5, yes = 1)) +
  scale_linetype_manual(values = c(no = 2, yes = 1)) +
  scale_x_continuous(limits = c(minplot, maxplot),
                     guide = guide_axis(check.overlap = TRUE)) +
  scale_y_continuous(limits = c(minplot, maxplot),
                     guide = guide_axis(check.overlap = TRUE)) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

tiff("C:/Users/samvida/Documents/Lindgren Group/Adiposity_Primary_Care/Reports/Manuscript/figures/abdominal_obesity_effects/WCadjBMI_vs_BMI.tiff",
     height = 3.5, width = 3.5, units = "cm",
     res = 300)
print(wcadjbmi_plot)
dev.off()


# Test sex heterogeneity ----

full_dat <- dat %>% 
  filter(pheno_tested %in% c("WCadjBMI", "WHRadjBMI") & sex_strata != "sex_comb") %>%
  select(all_of(c("beta", "se", "pval", "sex_strata", "pheno_tested", "SNP"))) 

full_dat <- split(full_dat, f = full_dat$pheno_tested)

het_results <- lapply(c("WCadjBMI", "WHRadjBMI"), function (p) {
  res <- full_dat[[p]] %>% pivot_wider(id_cols = SNP,
                                       names_from = sex_strata,
                                       values_from = c(beta, se, pval))
  
  res <- res %>%
    mutate(het_zstat = (beta_F - beta_M)/sqrt(se_F^2 + se_M^2),
           het_pval = pnorm(het_zstat, 0, 1, lower.tail = T),
           pheno_tested = p)
  return(res)
})
het_results <- bind_rows(het_results)

