# Author: Samvida S. Venkatesh
# Date: 08/06/2022

library(splines)
library(tidyverse)
theme_set(theme_bw())

plots_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/manuscript/figures/"

# For 3 colours that are sequential gray:
custom_grey_sequential <- c("#242424", "#7D7D7D", "#AEAEAE")

VARIDS <- c("rs429358", "rs8047587")

# Read data ----

# Genotypes / dosages at variants of interst
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

# Convert dosages to 0/1/2 genotypes based on threshold
DOSAGE_THRESHOLD_0 <- 0.5
DOSAGE_THRESHOLD_2 <- 1.5
var_dosages <- lapply(VARIDS, function (varid) {
  res <- var_dosages[[varid]] %>% 
    mutate(genotype = ifelse(!!as.symbol(varid) < DOSAGE_THRESHOLD_0, "0",
                             ifelse(!!as.symbol(varid) > DOSAGE_THRESHOLD_2, 
                                    "2", "1")),
           eid = as.character(eid)) %>%
    select(all_of(c("eid", "genotype")))
  return (res)
})
names(var_dosages) <- VARIDS

# Plot results ----

model_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/results/fit_objects_BMI_sex_comb.rds")

B <- model_dat$B
spline_posteriors <- model_dat$spline_posteriors
model_resid_var <- model_dat$resid_var

# Create mean coefficient matrix
# Matrix of means (id x basis)
mn_mat <- lapply(spline_posteriors, function (spobj) {
  return (as.data.frame(t(spobj$mu)))
})
mn_mat <- bind_rows(mn_mat)
rownames(mn_mat) <- names(spline_posteriors)
# Baseline data
mn_mat <- mn_mat - mn_mat[, 1]

# Get mean and S.E. of coefficients for ids by genotype 
calcMeanSDCoefs <- function (id_class) {
  sub_mn <- mn_mat
  sub_mn$eid <- as.character(rownames(sub_mn))
  sub_mn <- left_join(sub_mn, id_class, by = "eid")
  sub_mn <- sub_mn %>% filter(!is.na(sub_mn$genotype))
  
  # Calculate mean and s.e. of coefs
  coef_means <- sub_mn %>% 
    group_by(genotype) %>% summarise(across(-eid, mean))
  coef_sds <- sub_mn %>% 
    group_by(genotype) %>% summarise(across(-eid, sd))
  
  geno_counts <- sub_mn %>% group_by(genotype) %>% summarise(count = n())
  coef_ses <- coef_sds
  for (gi in unique(geno_counts$genotype)) {
    get_row <- which(coef_sds$genotype == gi)
    coef_ses[get_row, -1] <- coef_sds[get_row, -1] / sqrt(geno_counts$count[geno_counts$genotype == gi])
  }
  
  return (list(coef_means = coef_means, coef_ses = coef_ses))
}

## Predictions for a given coefficient matrix
getPredValues <- function (coef_mat, genos) {
  pred_vals <- t(apply(coef_mat, 1, 
                       function (x) B %*% x))
  # Wrangle into ggplot format
  for_plot <- as.data.frame(pred_vals)
  colnames(for_plot) <- paste0("d", 1:ncol(for_plot))
  
  for_plot$genotype <- factor(as.character(genos))
  
  for_plot <- for_plot %>% pivot_longer(cols = -genotype,
                                        names_to = "t_diff", 
                                        names_prefix = "d", 
                                        values_to = "pred_value") %>%
    mutate(t_diff = as.numeric(t_diff))
  return (for_plot)
}

## Cluster centroid mean and S.E., wrangled into long format for plot
getPlotDatGenotypes <- function (id_class) {
  
  # Get mean and s.d. coefficient matrices
  var_coefs <- calcMeanSDCoefs(id_class)
  genos <- var_coefs$coef_means[, 1]$genotype
  mn_coefs <- var_coefs$coef_means[, -1]
  lci_coefs <- var_coefs$coef_means[, -1] - 1.96*var_coefs$coef_ses[, -1]
  uci_coefs <- var_coefs$coef_means[, -1] + 1.96*var_coefs$coef_ses[, -1]
  
  # Predict full trajectory (returns centroid x t_diff matrix)
  pred_mns <- getPredValues(mn_coefs, genos) %>%
    rename(pred_mean_value = pred_value)
  pred_locis <- getPredValues(lci_coefs, genos) %>%
    rename(pred_loci_value = pred_value)
  pred_upcis <- getPredValues(uci_coefs, genos) %>%
    rename(pred_upci_value = pred_value)
  
  # Wrangle into ggplot format
  for_plot <- full_join(pred_mns, pred_locis, by = c("genotype", "t_diff"))
  for_plot <- full_join(for_plot, pred_upcis, by = c("genotype", "t_diff"))
  
  for_plot <- for_plot %>%
    mutate(t_diff_yrs = (t_diff - 1)/365)
  
  return (for_plot)
}

## Apply ----

lapply(VARIDS, function (v) {
    plot_dat <- getPlotDatGenotypes(var_dosages[[v]])
    
    res_plot <- ggplot(plot_dat, aes(x = t_diff_yrs, 
                                     y = pred_mean_value,
                                     col = genotype, fill = genotype)) +
      geom_line() +
      geom_ribbon(aes(ymin = pred_loci_value, ymax = pred_upci_value),
                  alpha = 0.25, linetype = 0) +
      scale_color_manual(values = custom_grey_sequential, guide = "none") +
      scale_fill_manual(values = custom_grey_sequential, guide = "none") +
      scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text = element_text(size = 8))
    
    tiff(filename = paste0(plots_dir, v, "_modelled_BMI_sex_comb.tiff"),
         height = 4.5, width = 4.5, units = "cm",
         res = 300)
    print(res_plot)
    dev.off()
})

