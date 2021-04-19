# Author: Samvida S. Venkatesh
# Date: 08/04/21

library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

# Read data ----

# Slopes
adj_slopes <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/adiposity_adj_slopes.rds")
# Raw data
adiposity <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/stratified_adiposity.rds")
# EID x disease matrix
eid_pheno_matrix <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/eid_phenotype_matrix.txt", 
                               sep = "\t", header = T, stringsAsFactors = F)
rownames(eid_pheno_matrix) <- eid_pheno_matrix$eid
eid_pheno_matrix <- eid_pheno_matrix[, -1]
# Disease dictionary
dictionary <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt",
                         sep = "\t", header = T, stringsAsFactors = F)
DISEASES <- dictionary$phenotype
colnames(eid_pheno_matrix) <- DISEASES

PHENOTYPES <- c("BMI", "WHR")
STRATA <- c("white_F", "white_M", "white_sexcomb")

# Calculate overlaps between BMI and WHR ----

two_by_two <- lapply(STRATA, function (s) {
  WHR_df <- adj_slopes[["WHR"]][[s]]
  BMI_df <- adj_slopes[["BMI"]][[s]]
  
  BMIg <- BMI_df$eid[BMI_df$gainer]
  BMIl <- BMI_df$eid[!BMI_df$gainer]
  
  WHRg <- WHR_df$eid[WHR_df$gainer]
  WHRl <- WHR_df$eid[!WHR_df$gainer]
  
  gg <- intersect(BMIg, WHRg)
  gl <- intersect(BMIg, WHRl)
  lg <- intersect(BMIl, WHRg)
  ll <- intersect(BMIl, WHRl)
  all_overlap <- c(gg, gl, lg, ll)
  
  # Annotate dfs with group belonging
  BMI_df <- subset(BMI_df, BMI_df$eid %in% all_overlap)
  WHR_df <- subset(WHR_df, WHR_df$eid %in% all_overlap)
  
  BMI_df$overlap_group <- ifelse(BMI_df$eid %in% gg, "BMIg-WHRg",
                                 ifelse(BMI_df$eid %in% gl, "BMIg-WHRl",
                                        ifelse(BMI_df$eid %in% lg, "BMIl-WHRg",
                                               ifelse(BMI_df$eid %in% ll, "BMIl-WHRl",
                                                      "error"))))
  
  WHR_df$overlap_group <- ifelse(WHR_df$eid %in% gg, "BMIg-WHRg",
                                 ifelse(WHR_df$eid %in% gl, "BMIg-WHRl",
                                        ifelse(WHR_df$eid %in% lg, "BMIl-WHRg",
                                               ifelse(WHR_df$eid %in% ll, "BMIl-WHRl",
                                                      "error"))))
  
  # Print percent of IDs in each group
  tot <- length(all_overlap)
  tab <- data.frame(BMIg = c(length(gg)/tot, length(gl)/tot), 
                    BMIl = c(length(lg)/tot, length(ll)/tot))
  rownames(tab) <- c("WHRg", "WHRl")
  print(tab)
  
  return (list(BMI = BMI_df, WHR = WHR_df))
})
names(two_by_two) <- STRATA

# Mean trajectories for individuals in overlap groups ----

plot_mean_traj <- lapply(STRATA, function (s) {
  p_dat <- lapply(PHENOTYPES, function (p) {
    df <- inner_join(adiposity[[p]][[s]], 
                     two_by_two[[s]][[p]][, c("eid", "overlap_group")], 
                     by = "eid")
    # Calculate mean and SE in each 5-year interval within each overlap gp
    age_bin_cuts <- seq(20, 80, by = 5)
    df$age_bin <- cut(df$age_event, age_bin_cuts, include.lowest = T)
    p_df <- df %>% group_by(overlap_group, age_bin) %>% 
      summarise(count = n(),
                mean_value = mean(value),
                se_value = sd(value)/sqrt(count))
    p_df$adipo_trait <- p
    return (p_df)
  })
  p_dat <- bind_rows(p_dat)
  # Plot 
  p <- ggplot(p_dat, aes(x = age_bin, y = mean_value,
                         group = overlap_group, 
                         color = overlap_group, 
                         fill = overlap_group)) +
    facet_wrap(~adipo_trait, nrow = 2, scales = "free_y") + 
    geom_point() +
    geom_path() +
    geom_ribbon(aes(ymin = mean_value - se_value, 
                    ymax = mean_value + se_value),
                alpha = 0.2) +
    scale_fill_manual(values = c("BMIg-WHRg" = "#984EA3",
                                 "BMIg-WHRl" = "#4DAF4A",
                                 "BMIl-WHRg" = "#377EB8",
                                 "BMIl-WHRl" = "#E41A1C")) +
    scale_color_manual(values = c("BMIg-WHRg" = "#984EA3",
                                  "BMIg-WHRl" = "#4DAF4A",
                                  "BMIl-WHRg" = "#377EB8",
                                  "BMIl-WHRl" = "#E41A1C")) +
    labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity",
         title = s)
  return (p)
})

pdf("plots/multivariate/BMI_WHR_groups_trajectories.pdf")
plot_mean_traj
dev.off()


# Disease enrichment in groups ---- 

## Calculate code counts and frequencies in each group and full cohort ----

counts <- lapply(STRATA, function (s) {
  
  df <- two_by_two[[s]][[1]]
  
  gg <- df$eid[df$overlap_group == "BMIg-WHRg"]
  gl <- df$eid[df$overlap_group == "BMIg-WHRl"]
  lg <- df$eid[df$overlap_group == "BMIl-WHRg"]
  ll <- df$eid[df$overlap_group == "BMIl-WHRl"]
  all_overlap <- c(gg, gl, lg, ll)
  
  cohort_counts <- colSums(eid_pheno_matrix[all_overlap, ], na.rm = T)
  gg_counts <- colSums(eid_pheno_matrix[gg, ], na.rm = T)
  gl_counts <- colSums(eid_pheno_matrix[gl, ], na.rm = T)
  lg_counts <- colSums(eid_pheno_matrix[lg, ], na.rm = T)
  ll_counts <- colSums(eid_pheno_matrix[ll, ], na.rm = T)
  
  counts <- data.frame(disease = DISEASES, cohort_counts, 
                       gg_counts, gl_counts, lg_counts, ll_counts,
                       NCOHORT = length(df$eid),
                       NGG = length(gg), NGL = length(gl),
                       NLG = length(lg), NLL = length(ll))
  
  # Convert counts to frequencies
  counts$cohort_freq <- counts$cohort_counts / counts$NCOHORT
  counts$gg_freq <- counts$gg_counts / counts$NGG
  counts$gl_freq <- counts$gl_counts / counts$NGL
  counts$lg_freq <- counts$lg_counts / counts$NLG
  counts$ll_freq <- counts$ll_counts / counts$NLL
  
  return (counts)
  
})

names(counts) <- STRATA

## Enrichment tests for codes in cases vs controls ----

enrichment_tests <- lapply(STRATA, function (s) {
  df <- counts[[s]]
  
  # Fisher's test for overrepresentation applied to each disease
  pval_gg <- rep(NA, nrow(df))
  pval_gl <- rep(NA, nrow(df))
  pval_lg <- rep(NA, nrow(df))
  pval_ll <- rep(NA, nrow(df))
  
  for (i in 1:nrow(df)) {
    r <- df[i, ]
    # GG compared to cohort
    mat_gg <- matrix(c(r[["gg_counts"]],
                       r[["cohort_counts"]],
                       df$NGG[1] - r[["gg_counts"]],
                       df$NCOHORT[1] - r[["cohort_counts"]]), 2, 2)
    pval_gg[i] <- fisher.test(mat_gg)$p.value
    # GL compared to cohort
    mat_gl <- matrix(c(r[["gl_counts"]],
                       r[["cohort_counts"]],
                       df$NGL[1] - r[["gl_counts"]],
                       df$NCOHORT[1] - r[["cohort_counts"]]), 2, 2)
    pval_gl[i] <- fisher.test(mat_gl)$p.value
    # LG compared to cohort
    mat_lg <- matrix(c(r[["lg_counts"]],
                       r[["cohort_counts"]],
                       df$NLG[1] - r[["lg_counts"]],
                       df$NCOHORT[1] - r[["cohort_counts"]]), 2, 2)
    pval_lg[i] <- fisher.test(mat_lg)$p.value
    # LL compared to cohort
    mat_ll <- matrix(c(r[["ll_counts"]],
                       r[["cohort_counts"]],
                       df$NLL[1] - r[["ll_counts"]],
                       df$NCOHORT[1] - r[["cohort_counts"]]), 2, 2)
    pval_ll[i] <- fisher.test(mat_ll)$p.value
  }
  
  df$pval_gg <- pval_gg
  df$pval_gl <- pval_gl
  df$pval_lg <- pval_lg
  df$pval_ll <- pval_ll
  
  df$FDR_gg <- p.adjust(df$pval_gg, method = "fdr")
  df$FDR_gl <- p.adjust(df$pval_gl, method = "fdr")
  df$FDR_lg <- p.adjust(df$pval_lg, method = "fdr")
  df$FDR_ll <- p.adjust(df$pval_ll, method = "fdr")
  
  df$effect_gg <- ifelse(df$gg_freq > df$cohort_freq,
                         "over-represented",
                         "under-represented")
  df$effect_gl <- ifelse(df$gl_freq > df$cohort_freq,
                         "over-represented",
                         "under-represented")
  df$effect_lg <- ifelse(df$lg_freq > df$cohort_freq,
                         "over-represented",
                         "under-represented")
  df$effect_ll <- ifelse(df$ll_freq > df$cohort_freq,
                         "over-represented",
                         "under-represented")
  
  write.table(df, paste0("results/phenotype_analysis/BMI_WHR_overlap_disease_enrichment_", 
                         s, ".txt"), 
              sep = "\t", quote = F, row.names = F)
  
  return (df)
})
names(enrichment_tests) <- STRATA

saveRDS(enrichment_tests, "/well/lindgren/UKBIOBANK/samvida/adiposity/results/phenotype_analysis/BMI_WHR_overlap_disease_enrichment.rds")
