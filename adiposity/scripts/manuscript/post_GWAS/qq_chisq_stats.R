# Author: Samvida S. Venkatesh
# Date: 12/04/2023

library(tidyverse)
theme_set(theme_bw())

# Read data ----

mainpath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/enrichment_change_signal_baseline"
PHENO <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")
PARAMETERS <- c("b1", "k1", "k1_k2", "k1_k2_k3")

sumstats <- lapply(PHENO, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    per_par <- lapply(PARAMETERS, function (pr) {
      res <- read.table(paste0(mainpath, "/", p, "_", sx, "_", pr, "_grepped_sumstats.txt"), 
                        sep = "\t", header = T, 
                        comment.char = "@", stringsAsFactors = F)
      res$pheno <- p
      res$sex_strata <- sx
      res$parameter <- pr
      return (res)
    })
    names(per_par) <- PARAMETERS
    return (per_par)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(sumstats) <- PHENO

# Plot QQ-plots of observed vs expected chi-square statistics ----

getExpChiSq <- function (obs_chisq) {
  obs_chisq <- sort(obs_chisq)
  # Get expected chi-squared stats
  k <- 1:length(obs_chisq)
  beta_k <- k/(length(obs_chisq)+1)
  exp_chisq <- qchisq(beta_k, df = 1)
  
  # CI based on the beta distribution
  # parameters: a = k, b = n-k where k=seq(1, n)
  beta_lci <- qbeta(0.025, k, rev(k))
  beta_uci <- qbeta(0.975, k, rev(k))
  # Generate the chi-square quantile from these probabilities
  exp_lci <- qchisq(beta_lci, df = 1)
  exp_uci <- qchisq(beta_uci, df = 1)
  
  plot_qq <- data.frame(obs = obs_chisq, exp = exp_chisq,
                        lci = exp_lci, uci = exp_uci)
  return (plot_qq)
}

plotQQ <- function (qq_dat) {
  plot_res <- ggplot(qq_dat, 
                     aes(x = exp, y = obs)) +
    geom_ribbon(aes(ymin = lci, ymax = uci), 
                fill = "grey", alpha = 0.5) +
    geom_point(size = 0.7) +
    geom_abline(intercept = 0, slope = 1) +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 6))
  return (plot_res)
}

lapply(PHENO, function (p) {
  lapply(SEX_STRATA, function (sx) {
    lapply(PARAMETERS, function (pr) {
      qq_dat <- getExpChiSq(obs_chisq = sumstats[[p]][[sx]][[pr]]$CHISQ_BOLT_LMM_INF)
      
      ggsave(paste0(mainpath, "/plots/", p, "_", sx, "_", pr, "_qq_chisq.png"),
             plotQQ(qq_dat),
             units = "cm", height = 5, width = 5)
    })
  })
})

# Plot histogram of Z-score(b0)*Z-score(b1) ----


