# Author: Samvida S. Venkatesh
# Date: 29/06/21

library(tidyverse)
theme_set(theme_bw())
library(RColorBrewer)

PHENO <- read.table("pheno_names.txt", header = F)$V1
SEX_STRATA <- c("F", "M", "sex_comb")

# Read files ----

# EID x disease matrix
eid_pheno_matrix <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/eid_phenotype_matrix.txt", 
                               sep = "\t", header = T, stringsAsFactors = F)
rownames(eid_pheno_matrix) <- eid_pheno_matrix$eid
eid_pheno_matrix <- eid_pheno_matrix[, -1]

# Disease dictionary
dictionary <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt",
                         sep = "\t", header = T, stringsAsFactors = F, quote = "")
DISEASES <- dictionary$phenotype
colnames(eid_pheno_matrix) <- DISEASES

# Cluster assignment
cdat <- lapply(PHENO, function (pheno) {
  res <- lapply(SEX_STRATA, function (sx) {
    df <- read.table(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/splines/new_models/",
                      pheno, "/coefficients_", pheno, "_", sx, ".txt"),
               header = T, sep = "\t", stringsAsFactors = F)
    return (df)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(cdat) <- PHENO

# Covariates
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/model_covariates.rds")
PCs <- paste0("PC", 1:21)

# Logistic regression functions ----

prep_indiv_data <- function (pheno, sx) {
  # Get cluster information
  dat_for_test <- cdat[[pheno]][[sx]]
  # Code cluster variables as dummy variables for regression
  dat_for_test$dummy <- 1
  dat_for_test <- pivot_wider(dat_for_test, id_cols = eid, 
                              names_from = cluster, names_prefix = "c",
                              values_from = dummy, values_fill = 0)
  # Add in covariate information
  dat_for_test <- merge(dat_for_test, covars[[pheno]], by = "eid")
  return (dat_for_test)
}

get_formula <- function (nclust, sx, adj_baseline = T) {
  clust_vars <- paste0("c", 2:nclust)
  test_form <- paste0("disease_status ~ ", 
                      paste(PCs, collapse = " + "),
                      " + ",
                      paste(clust_vars, collapse = " + "))
  # Add covariate for sex in sex-combined analysis
  if (sx == "sex_comb") {
    test_form <- paste0(test_form, " + sex")
  }
  # Add covariate for baseline trait in adj- analysis
  if (adj_baseline) {
    test_form <- paste0(test_form, " + baseline_trait")
  }
  return (test_form)
}

prep_dis_data <- function (ids) {
  # Subset disease matrix to relevant subset
  dis_mat <- eid_pheno_matrix[as.character(ids), ]
  # Only keep diseases that have at least one patient
  # represented in the cohort
  keep_diseases <- colSums(dis_mat, na.rm = T) > 0 & 
    colSums(dis_mat, na.rm = T) < nrow(dis_mat)
  dis_mat <- dis_mat[, keep_diseases]
  return (dis_mat)
}

perform_logreg <- function (pheno, sx, adj_baseline) {
  
  indiv_dat <- prep_indiv_data(pheno, sx)
  NCLUST <- length(unique(cdat[[pheno]][[sx]]$cluster))
  clust_covars <- paste0("c", 2:NCLUST)

  test_form <- get_formula(NCLUST, sx, adj_baseline)
  
  dis_mat <- prep_dis_data(indiv_dat$eid)
  DISEASES <- colnames(dis_mat)
  dis_mat$eid <- rownames(dis_mat)
  
  # Loop through each disease and apply formula
  res <- lapply(DISEASES, function (d) {
    test_dis <- dis_mat[, c("eid", d)]
    colnames(test_dis)[2] <- "disease_status"
    test_df <- merge(indiv_dat, test_dis, by = "eid")
    test_res <- summary(glm(as.formula(test_form), family = binomial,
                   data = test_df))
    res_df <- data.frame(cluster = clust_covars,
                         beta = test_res$coefficients[clust_covars, 1],
                         se = test_res$coefficients[clust_covars, 2],
                         pval = test_res$coefficients[clust_covars, 4],
                         disease = d)
    return (res_df)
  })
  # Bind all disease results and retain only terms that 
  # pass FDR adjustment (FDR < 0.05)
  res <- bind_rows(res)
  res$FDR <- p.adjust(res$pval, "fdr")
  return (res)
}

# Apply tests across all phenotypes and sexes ----

logreg_tests <- lapply(PHENO, function (pheno) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    # Unadjusted for baseline trait value
    unadj <- perform_logreg(pheno, sx, adj_baseline = F)
    unadj$baseline_trait_adj <- F
    # Adjustment for baseline trait value
    adj <- perform_logreg(pheno, sx, adj_baseline = T)
    adj$baseline_trait_adj <- T
    res <- bind_rows(unadj, adj)
    write.table(res, paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/splines/new_models/",
                              pheno, "/cluster_disease_logreg_", 
                              pheno, "_", sx, ".txt"),
                sep = "\t", quote = F, row.names = F)
    return (res)
  })
  names(res_list) <- SEX_STRATA
})
names(logreg_tests) <- PHENO
saveRDS(logreg_tests, 
        "/well/lindgren/UKBIOBANK/samvida/adiposity/splines/new_models/all_phenotypes_disease_logreg.rds")

# Cohort and cluster disease prevalence ----

prevalence_tables <- lapply(PHENO, function (pheno) {
  res_tables <- lapply(SEX_STRATA, function (sx) {
    # Get relevant disease matrix
    relev_cdat <- cdat[[pheno]][[sx]]
    summ_mat <- eid_pheno_matrix[as.character(relev_cdat$eid), ]
    # Get cohort frequencies
    cohort_count <- nrow(relev_cdat)
    cohort_freqs <- colSums(summ_mat) / cohort_count
    # Group by cluster and get cluster frequencies
    summ_mat$eid <- rownames(summ_mat)
    summ_mat$cluster <- relev_cdat$cluster[match(relev_cdat$eid,
                                                 summ_mat$eid)]
    clust_counts <- relev_cdat %>% count(cluster)
    summ_mat <- summ_mat %>% group_by(cluster) %>%
      summarise(across(all_of(DISEASES), sum))
    summ_mat <- t(summ_mat)
    summ_mat <- as.data.frame(t(apply(summ_mat[DISEASES, ], 1, 
                        function (x) x / clust_counts$n)))
    colnames(summ_mat) <- paste0("cluster_", 1:ncol(summ_mat))
    summ_mat$cohort_freq <- cohort_freqs
    return (summ_mat)
  })
  names(res_tables) <- SEX_STRATA
  return (res_tables)
})
names(prevalence_tables) <- PHENO

saveRDS(prevalence_tables, "all_pheno_cluster_disease_prevalence_tables.rds")

# Fisher's enrichment test ----

# Returns table of case and control frequencies, p-values for each disease

enrich_test <- function (control_ids, case_ids) {
  
  # Subset disease matrix to relevant results ----
  dis_mat <- eid_pheno_matrix[rownames(eid_pheno_matrix) %in% 
                                c(control_ids, case_ids), ]
  # Only keep diseases that have at least one patient
  # represented in the cohort
  keep_diseases <- colSums(dis_mat, na.rm = T) > 0 & 
    colSums(dis_mat, na.rm = T) < nrow(dis_mat)
  dis_mat <- dis_mat[, keep_diseases]
  DISEASES <- colnames(dis_mat)
  
  # Get disease counts and frequency in controls and cases ----
  control_counts <- colSums(dis_mat[control_ids, ], na.rm = T)
  case_counts <- colSums(dis_mat[case_ids, ], na.rm = T)
  
  counts <- data.frame(disease = DISEASES,
                       control_counts, case_counts,
                       NCONTROL = length(control_ids),
                       NCASES = length(case_ids))
  counts$control_freq <- counts$control_counts / counts$NCONTROL
  counts$case_freq <- counts$case_counts / counts$NCASES
  
  # Test for enrichment ----
  pvals <- rep(NA, length(DISEASES))
  for (d in 1:length(DISEASES)) {
    r <- counts[d, ]
    test_mat <- matrix(c(r[["case_counts"]], 
                         r[["control_counts"]],
                         r[["NCASES"]] - r[["case_counts"]],
                         r[["NCONTROL"]] - r[["control_counts"]]), 2, 2)
    pvals[d] <- fisher.test(test_mat)$p.value
  }
  
  # Report results ----
  res <- counts
  res$pval <- pvals
  res$effect_case <- ifelse(res$case_freq > res$control_freq, 
                            "over-represented", 
                            ifelse(res$case_freq < res$control_freq,
                                   "under-represented", "no change"))
  res$FDR <- p.adjust(res$pval)
  res <- res[, c("disease", "NCASES", "case_counts",
                 "pval", "FDR", "effect_case")]
  return (res)
}

# Apply enrichment test to each cluster in each strata ----

disease_tables <- lapply(PHENO, function (pheno) {
  res <- lapply(SEX_STRATA, function (sx) {
    id_dat <- cdat[[pheno]][[sx]]
    NCLUST <- length(unique(id_dat$cluster))
    
    per_clust_res <- lapply(1:NCLUST, function (k) {
      k_ids <- id_dat$eid[id_dat$cluster == k]
      k_ids <- as.character(k_ids)
      control_ids <- id_dat$eid[id_dat$cluster != k]
      control_ids <- as.character(control_ids)
      test_res <- enrich_test(control_ids = control_ids, 
                              case_ids = k_ids)
      colnames(test_res)[-1] <- paste0("cluster", k, "_", 
                                       colnames(test_res)[-1])
      return (test_res)
    })
    all_res <- per_clust_res %>% reduce(full_join, by = "disease")
    # For results, report size of cohort and disease frequency in cohort
    all_res$NCOHORT <- nrow(id_dat)
    # Get disease counts
    dcohort <- colSums(eid_pheno_matrix[as.character(id_dat$eid), 
                                        all_res$disease],
                       na.rm = T)
    all_res$cohort_freq <- dcohort
    write.table(all_res, paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/splines/new_models/",
                       pheno, "/cluster_disease_enrichment_", 
                       pheno, "_", sx, ".txt"),
                sep = "\t", quote = F, row.names = F)
    return (all_res)
  })
  names(res) <- SEX_STRATA
})
names(disease_tables) <- PHENO

saveRDS(disease_tables, "all_phenotypes_disease_tables.rds")


