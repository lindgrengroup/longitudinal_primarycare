# Author: Samvida S. Venkatesh
# Date: 31/08/22

library(tidyverse)
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

vars_to_replicate <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/data/adipo_change_snps_replicate.txt",
                                sep = "\t", header = T, stringsAsFactors = F)
VARIDS <- vars_to_replicate$SNP

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
                paste0("PC", 1:21))
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

quantTest <- function (dat, ss) {
  covars_include <- MOD_COVARS
  if (ss != "sex_comb") {
    dat <- dat %>% filter(sex == ss)
    covars_include <- covars_include[-which(covars_include == "sex")]
  }
  if (length(unique(dat$data_provider)) < 2)
    covars_include <- covars_include[-which(covars_include == "data_provider")]
  
  mod_formula <- paste0("value_diff ~ dosage + ",
                        paste0(covars_include, collapse = " + "))
  
  modeled_dat <- lm(formula(mod_formula), data = dat)
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
                      "year_of_birth", "data_provider", 
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
      labs(x = paste0("Change in ", qp), 
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
    
    # Run regressions and plots
    regression_res <- lapply(SEX_STRATA, function (ss) {
      per_visit <- lapply(c(2, 3), function (vc) {
        sub_dat <- model_dat %>% filter(visit == vc)
        res <- quantTest(sub_dat, ss) %>%
          mutate(visit_compared = vc)
      })
      per_visit <- bind_rows(per_visit) %>%
        mutate(sex_strata = ss)
      
      return (per_visit)
    })
    regression_res <- bind_rows(regression_res) %>%
      mutate(pheno_tested = qp)
    
    
    # Plot results
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

