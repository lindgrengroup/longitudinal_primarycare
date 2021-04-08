# Author: Samvida S. Venkatesh
# Date: 26/03/2021

library(tidyverse)
theme_set(theme_bw())
library(RColorBrewer)
library(ggrepel)

# Read data ----

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

# Final slopes
final_slopes <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/adiposity_adj_slopes.rds")
PHENOTYPES <- names(final_slopes)
STRATA <- names(final_slopes[[1]])

# Calculate code counts and frequencies in each group and full cohort ----

counts <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    df <- final_slopes[[p]][[s]]
    df$q <- cut(df$residual, quantile(df$residual), include.lowest = T,
                labels = paste0("q", 1:4))
    
    gs <- df$eid[df$gainer]
    q1 <- df$eid[df$q == "q1"]
    q2 <- df$eid[df$q == "q2"]
    q3 <- df$eid[df$q == "q3"]
    q4 <- df$eid[df$q == "q4"]
    
    cohort_counts <- colSums(eid_pheno_matrix[df$eid, ], na.rm = T)
    gainer_counts <- colSums(eid_pheno_matrix[gs, ], na.rm = T)
    q1_counts <- colSums(eid_pheno_matrix[q1, ], na.rm = T)
    q2_counts <- colSums(eid_pheno_matrix[q2, ], na.rm = T)
    q3_counts <- colSums(eid_pheno_matrix[q3, ], na.rm = T)
    q4_counts <- colSums(eid_pheno_matrix[q4, ], na.rm = T)
    
    counts <- data.frame(disease = DISEASES,
                         cohort_counts, gainer_counts,
                         q1_counts, q2_counts,
                         q3_counts, q4_counts,
                         NCOHORT = length(df$eid),
                         NGAINERS = length(gs),
                         NQ1 = length(q1), NQ2 = length(q2),
                         NQ3 = length(q3), NQ4 = length(q4))
    
    # Convert counts to frequencies
    counts$cohort_freq <- counts$cohort_counts / counts$NCOHORT
    counts$gainer_freq <- counts$gainer_counts / counts$NGAINERS
    counts$q1_freq <- counts$q1_counts / counts$NQ1
    counts$q2_freq <- counts$q2_counts / counts$NQ2
    counts$q3_freq <- counts$q3_counts / counts$NQ3
    counts$q4_freq <- counts$q4_counts / counts$NQ4
    
    return (counts)
    
  })
  names(res) <- STRATA
  return (res)
})
names(counts) <- PHENOTYPES

# Enrichment tests for codes in cases vs controls ----

enrichment_tests <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    df <- counts[[p]][[s]]
    
    # Fisher's test for overrepresentation applied to each disease
    pval_gainer <- rep(NA, nrow(df))
    pval_q1 <- rep(NA, nrow(df))
    pval_q2 <- rep(NA, nrow(df))
    pval_q3 <- rep(NA, nrow(df))
    pval_q4 <- rep(NA, nrow(df))
    
    for (i in 1:nrow(df)) {
      r <- df[i, ]
      # Gainers compared to cohort
      mat_gainer <- matrix(c(r[["gainer_counts"]], 
                             r[["cohort_counts"]],
                             df$NGAINERS[1] - r[["gainer_counts"]],
                             df$NCOHORT[1] - r[["cohort_counts"]]), 2, 2)
      pval_gainer[i] <- fisher.test(mat_gainer)$p.value
      # q1 compared to cohort
      mat_q1 <- matrix(c(r[["q1_counts"]],
                         r[["cohort_counts"]],
                         df$NQ1[1] - r[["q1_counts"]],
                         df$NCOHORT[1] - r[["cohort_counts"]]), 2, 2)
      pval_q1[i] <- fisher.test(mat_q1)$p.value
      # q2 compared to cohort
      mat_q2 <- matrix(c(r[["q2_counts"]],
                         r[["cohort_counts"]],
                         df$NQ2[1] - r[["q2_counts"]],
                         df$NCOHORT[1] - r[["cohort_counts"]]), 2, 2)
      pval_q2[i] <- fisher.test(mat_q2)$p.value
      # Q3 compared to cohort
      mat_q3 <- matrix(c(r[["q3_counts"]],
                         r[["cohort_counts"]],
                         df$NQ3[1] - r[["q3_counts"]],
                         df$NCOHORT[1] - r[["cohort_counts"]]), 2, 2)
      pval_q3[i] <- fisher.test(mat_q3)$p.value
      # q4 compared to cohort
      mat_q4 <- matrix(c(r[["q4_counts"]],
                         r[["cohort_counts"]],
                         df$NQ4[1] - r[["q4_counts"]],
                         df$NCOHORT[1] - r[["cohort_counts"]]), 2, 2)
      pval_q4[i] <- fisher.test(mat_q4)$p.value
    }
    df$pval_gainer <- pval_gainer
    df$pval_q1 <- pval_q1
    df$pval_q2 <- pval_q2
    df$pval_q3 <- pval_q3
    df$pval_q4 <- pval_q4
    
    df$FDR_gainer <- p.adjust(df$pval_gainer, method = "fdr")
    df$FDR_q1 <- p.adjust(df$pval_q1, method = "fdr")
    df$FDR_q2 <- p.adjust(df$pval_q2, method = "fdr")
    df$FDR_q3 <- p.adjust(df$pval_q3, method = "fdr")
    df$FDR_q4 <- p.adjust(df$pval_q4, method = "fdr")
    
    df$effect_gainer <- ifelse(df$gainer_freq > df$cohort_freq, 
                               "over-represented", 
                               "under-represented")
    df$effect_q1 <- ifelse(df$q1_freq > df$cohort_freq,
                           "over-represented",
                           "under-represented")
    df$effect_q2 <- ifelse(df$q2_freq > df$cohort_freq,
                           "over-represented",
                           "under-represented")
    df$effect_q3 <- ifelse(df$q3_freq > df$cohort_freq,
                           "over-represented",
                           "under-represented")
    df$effect_q4 <- ifelse(df$q4_freq > df$cohort_freq,
                           "over-represented",
                           "under-represented")
    return (df)
  })
  names(res) <- STRATA
  return (res)
})
names(enrichment_tests) <- PHENOTYPES

saveRDS(enrichment_tests, "/well/lindgren/UKBIOBANK/samvida/adiposity/results/phenotype_analysis/disease_enrichment.rds")

# Write table of enrichment test results ----

# Only for white-ancestry gainers in BMI and WHR
sub_pheno <- c("BMI", "WHR")
white_gainers <- lapply(sub_pheno, function (p) {
  white_strata <- STRATA[grep("white", STRATA)]
  df <- lapply(white_strata, function (s) {
    res <- enrichment_tests[[p]][[s]][, c("disease", 
                                          "cohort_counts", "gainer_counts",
                                          "NCOHORT", "NGAINERS", 
                                          "cohort_freq", "gainer_freq",
                                          "pval_gainer", "FDR_gainer", 
                                          "effect_gainer")]
    res$strata <- s
    res$adiposity_trait <- p
    write.table(res, 
                paste0("results/phenotype_analysis/disease_enrichment_", p, 
                       "_gainers_", s, ".txt"),
                sep = "\t", quote = F, row.names = F)
    return (res)
  })
  names(df) <- white_strata
  # Summary of disease-wise results across strata
  df_summ <- bind_rows(df)
  df_summ <- pivot_wider(df_summ, id_cols = disease,
                         values_from = c(cohort_counts, gainer_counts,
                                         NCOHORT, NGAINERS, 
                                         cohort_freq, gainer_freq,
                                         pval_gainer, FDR_gainer, effect_gainer),
                         names_from = strata, names_sep = "_")
  
  write.table(df_summ, 
              paste0("results/phenotype_analysis/disease_enrichment_summary_", p, 
                     "_gainers.txt"),
              sep = "\t", quote = F, row.names = F)
  return (df)
})
names(white_gainers) <- sub_pheno

# Manhattan-style plots ----

# To mirror axis, write a rough script to move the x-axis to y = 0
# Adapted from: https://stackoverflow.com/questions/39071002/moving-x-or-y-axis-together-with-tick-labels-to-the-middle-of-a-single-ggplot-n

shift_x_axis <- function(p, y = 0) {
  g <- ggplotGrob(p)
  dummy <- data.frame(y = y)
  ax <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  p + annotation_custom(grid::grobTree(ax, 
                                       vp = grid::viewport(y=1, 
                                                           height=sum(ax$height))), 
                        ymax=y, ymin=y) +
    geom_hline(aes(yintercept=y), data = dummy) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
}

manhattan_plots <- lapply(PHENOTYPES, function (p) {
  
  res <- lapply(STRATA, function (s) {
    
    df <- enrichment_tests[[p]][[s]]
    df$plot_pval <- -log10(df$pval_gainer)
    
    # Flip "pval" for under-represented so it mirrors the y-axis
    df$plot_pval <- ifelse(df$effect_gainer == "under-represented",
                           -df$plot_pval, df$plot_pval)
    # Plotting parameters
    aticks <- round(seq(min(df$plot_pval), max(df$plot_pval),
                        length.out = 5), 2)
    alabels <- abs(aticks)
    # Name points with lowest p-values (or FDR < 0.01 if threshold too high)
    name_points_threshold <- max(quantile(abs(df$plot_pval), 0.95), 2)
    
    # Build plot
    p <- ggplot(df, aes(x = gainer_freq, y = plot_pval)) +
      geom_point(aes(col = effect_gainer)) +
      geom_text_repel(data = subset(df, abs(df$plot_pval) > name_points_threshold), 
                      aes(label = disease)) +
      scale_color_manual(values = c("#A50026", "#313695"),
                         guide = F) +
      scale_y_continuous(breaks = aticks, labels = alabels) +
      labs(y = "-log10(pval)", x = "Disease frequency in gainers",
           title = s)
    # Shift axes by custom function to mirror under-represented points
    p <- shift_x_axis(p, y = 0)
    
    return (p)
  })
  names(res) <- STRATA
  # Print plots in each trait
  pdf(paste0("plots/phenotype_analysis/enriched_diseases_gainers_", p, ".pdf"), 
      onefile = T)
  print(res)
  dev.off()
  return (res)
})
names(manhattan_plots) <- PHENOTYPES

