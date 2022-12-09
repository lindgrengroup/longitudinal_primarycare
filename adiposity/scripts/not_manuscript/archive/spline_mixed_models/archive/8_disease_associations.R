# Author: Samvida S. Venkatesh
# Date: 29/06/21

library(tidyverse)
theme_set(theme_bw())
library(RColorBrewer)

PHENO <- "BMI"
SEX_STRATA <- "F"

# Read files ----

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
slopes <- read.table(paste0("cluster_assignment_", PHENO, "_", SEX_STRATA, ".txt"),
           sep = "\t", header = T)
# The number of covariates may change so check this each time
RF_TERMS <- colnames(slopes)[32:(ncol(slopes) - 1)]

# Subset ID and disease matrix to relevant results ----- 

eid_pheno_matrix <- 
  eid_pheno_matrix[rownames(eid_pheno_matrix) %in% slopes$eid, ]

# Only keep diseases that have at least one case 
keep_diseases <- colSums(eid_pheno_matrix, na.rm = T) > 0 & 
  colSums(eid_pheno_matrix, na.rm = T) < nrow(eid_pheno_matrix)
eid_pheno_matrix <- eid_pheno_matrix[, keep_diseases]

DISEASES <- colnames(eid_pheno_matrix)

# Logistic regression for each of the random effect terms ----

# Prepare eidxpheno matrix
logit_mat <- eid_pheno_matrix
logit_mat$eid <- rownames(eid_pheno_matrix)

run_logit <- function (rf) {
  all_diseases <- lapply(DISEASES, function (dis) {
    tmp_df <- data.frame(eid = logit_mat$eid)
    tmp_df$disease_outcome <- factor(logit_mat[, dis])
    tmp_df$predictor <- slopes[match(slopes$eid, tmp_df$eid), rf]
    mod <- glm(disease_outcome ~ predictor, family = "binomial",
               data = tmp_df)
    res <- summary(mod)$coefficients["predictor", ]
    return (res)
  })
  all_diseases <- bind_rows(all_diseases)
  colnames(all_diseases) <- c("beta", "se", "z_stat", "p_value")
  # Convert beta and S.E. to OR and 95% CI
  all_diseases <- all_diseases %>% 
    mutate(OR = exp(beta), 
           lci = exp(beta - (1.96*se)),
           uci = exp(beta + (1.96*se)))
  all_diseases$DISEASE <- DISEASES
  # Write results table
  write.table(all_diseases, paste0("logit_regression_", rf, ".txt"),
              sep = "\t", row.names = F, quote = F)
  return (all_diseases)
}

rf_associations <- lapply(RF_TERMS, function (rf) run_logit (rf) )
names(rf_associations) <- RF_TERMS

## Plot associations ----

or_plots <- lapply(RF_TERMS, function (rf) {
  # Only plot top 25 diseases (in terms of p-value)
  plot_dat <- rf_associations[[rf]]
  plot_dat <- plot_dat[order(plot_dat$p_value), ]
  plot_dat <- plot_dat[1:min(25, nrow(plot_dat)), ]
  
  # Add marker for significance
  plot_dat$sig <- factor(ifelse(plot_dat$p_value < 0.05, 
                                "unadj. P < 0.05", "n.s."),
                         levels = c("unadj. P < 0.05", "n.s."))
  
  # Order by OR 
  plot_dat$DISEASE <- factor(plot_dat$DISEASE,
                             levels = plot_dat$DISEASE[order(plot_dat$OR)])
  
  res <- ggplot(plot_dat, aes(x = DISEASE, 
                              y = OR, ymin = lci, ymax = uci)) +
    geom_point(aes(alpha = sig), size = 2) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_errorbar(aes(ymin = lci, ymax = uci, linetype = sig,
                      alpha = sig)) +
    scale_y_log10(breaks = scales::pretty_breaks(n = 7)) +
    scale_x_discrete(label = function(x) stringr::str_trunc(x, 25)) +
    scale_alpha_discrete(range = c(1, 0.5)) +
    labs(title = paste0("Assoc. with ", rf)) +
    coord_flip() 
  return (res)
})

# Logistic regression for differences in the clusters ----

cluster_diffs <- lapply(DISEASES, function (dis) {
  tmp_df <- data.frame(eid = logit_mat$eid)
  tmp_df$disease_outcome <- factor(logit_mat[, dis])
  tmp_df$cluster <- factor(slopes$cluster[match(slopes$eid, tmp_df$eid)])
  mod <- glm(disease_outcome ~ cluster, family = "binomial",
             data = tmp_df)
  res <- data.frame(summary(mod)$coefficients)
  colnames(res) <- c("beta", "se", "z_stat", "p_value")
  res$cluster <- rownames(res)
  res$DISEASE <- dis
  res <- subset(res, res$p_value < 0.05)
  return (res)
})

cluster_diffs <- bind_rows(cluster_diffs)
cluster_diffs <- subset(cluster_diffs, 
                        cluster_diffs$cluster != "(Intercept)")

# Convert beta and S.E. to OR and 95% CI
cluster_diffs <- cluster_diffs %>% 
  mutate(OR = exp(beta), 
         lci = exp(beta - (1.96*se)),
         uci = exp(beta + (1.96*se)))

write.table(cluster_diffs, paste0("cluster_logit_regression.txt"),
            sep = "\t", row.names = F, quote = F)

## Plot association with clusters ----

# Only plot top 25 diseases (in terms of p-value)
plot_dat <- cluster_diffs
plot_dat <- plot_dat[order(plot_dat$p_value), ]
plot_dat <- plot_dat[1:min(25, nrow(plot_dat)), ]

# Order by OR 
plot_dat$DISEASE <- factor(plot_dat$DISEASE,
                           levels = plot_dat$DISEASE[order(plot_dat$OR)])

cluster_plot <- ggplot(plot_dat, aes(x = DISEASE, 
                            y = OR, ymin = lci, ymax = uci,
                            colour = cluster)) +
  facet_wrap(~cluster) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_errorbar(aes(ymin = lci, ymax = uci)) +
  scale_y_log10(breaks = scales::pretty_breaks(n = 7)) +
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 25)) +
  scale_colour_brewer(palette = "Dark2", guide = F) +
  labs(title = "Association with clusters") +
  coord_flip() 

pdf("logit_regression_plots.pdf", onefile = T)
print(or_plots)
print(cluster_plot)
dev.off()
