source("scripts/00_load_data.R")
list_import <- readRDS(file = paste0("output/linear_model_fitting/mat_vec_combine_task04b.RDS"))

control$N_LOCI <- 113
mu_beta <- rep(0, N_AGE_BIN)
n_param <- 2 + 2 * control$N_LOCI
N_AGE_BIN * n_param
str(list_import$mat)

Sig_beta <- Sig_eps_fn(len = 35, sig_GP = .5, sig_RI = .5, sig_resid = 1e-10)
# if (control$do_plots) {
beta_gen <- mvtnorm::rmvnorm(n = 10, mean = mu_beta, sigma = Sig_beta)
pdf("plots/test_beta_prior.pdf")
matplot(t(beta_gen), ty = "l")
dev.off()
# }


beta_post_var <- solve(kronecker(X = diag(rep(1, n_param)), Y = solve(Sig_beta)) + list_import$mat)
beta_post_mn <- beta_post_var %*% list_import$vec
beta_mat_mn<-beta_mat_sd <- matrix(NA, N_AGE_BIN, n_param)
beta_mat_mn[] <- c(beta_post_mn)
beta_mat_sd[] <- sqrt(diag(beta_post_var) * 1)
  
dir_gen <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/biobank_completion/BMI/variant_selection/bed_files"
gen_meta <- read.table(file.path(dir_gen, "BMI_filtered.lead_snps.txt"), header = TRUE)
rs_in_meta <- names(meta)[grepl("rs", names(meta))]
rs_fitted <- rs_in_meta[1:control$N_LOCI]

pdf("plots/beta_out.pdf", 12, 8)
par(mfrow = c(3, 4), oma = c(3, 3, 1, 1))
cextext <- .75
for (snpnum in 1:control$N_LOCI) {
  snpinds <- 2 + c(snpnum, snpnum + control$N_LOCI)
  gen_meta_curr <- gen_meta[match(rs_fitted[snpnum], gen_meta$rsid), ]
  beta_bmi_gwas <- gen_meta_curr$beta
  allele_freq <- 1 - gen_meta_curr$A1FREQ
  signed_beta_mat_mn <- beta_mat_mn#sign(beta_bmi_gwas) * beta_mat_mn * -1
  matplot(x = age_all,
          y = signed_beta_mat_mn[, snpinds], ty = "l", lty = 1, lwd = 3, 
          col = 1:2, ylab = "", xlab = "Age", las = 2,
          ylim = c(-1, 1) * max(abs(c(signed_beta_mat_mn[, snpinds] + 2 * beta_mat_sd[, snpinds],
                                  signed_beta_mat_mn[, snpinds] - 2 * beta_mat_sd[, snpinds]))))
  mtext(side = 3, text = paste0("freq ", gen_meta_curr$allele2, " = ", round(allele_freq, 2)), line = 0, cex = cextext)
  mtext(side = 3, text = paste0("beta = ", round(-beta_bmi_gwas, 3)), line = 1, cex = cextext)
  mtext(side = 3, text = rs_fitted[snpnum], line = 2, cex = cextext)
  mtext(side = 2, text = "log BMI", cex = cextext, line = 4)
  for (j in 1:2) {
    ind <- snpinds[j]
    lines(x = age_all, y = signed_beta_mat_mn[, ind] + 2 * beta_mat_sd[, ind], lty = 1, col = j)
    lines(x = age_all, y = signed_beta_mat_mn[, ind] - 2 * beta_mat_sd[, ind], lty = 1, col = j)
  }
  abline(h = 0)
  abline(h = seq(-.1, .1, by = .01), col = 1, lty = 3)
  legend(x = ifelse(-beta_bmi_gwas < 0, "topleft", "bottomleft"), legend = c(1, 2), 
         title = paste("# copies ", gen_meta_curr$allele2), col = 1:2, lty = 1, bg = "white")
}
dev.off()




pdf("plots/sex_beta_out.pdf", 12, 8)
par(mfrow = c(1, 1))
snpinds <- 1:2
matplot(x = age_all,
        y = beta_mat_mn[, snpinds], ty = "l", lty = 1, lwd = 3, col = 1:2, ylab = "log BMI", xlab = "Age",
         ylim = range((c(beta_mat_mn[, snpinds] + 2 * beta_mat_sd[, snpinds],
                                          beta_mat_mn[, snpinds] - 2 * beta_mat_sd[, snpinds]))))
for (j in 1:2) {
  ind <- snpinds[j]
  lines(x = age_all, y = beta_mat_mn[, ind] + 2 * beta_mat_sd[, ind], lty = 1, col = j)
  lines(x = age_all, y = beta_mat_mn[, ind] - 2 * beta_mat_sd[, ind], lty = 1, col = j)
}
legend(x = "bottomright", legend = c("Female", "Male"), col = 1:2, lty = 1)
abline(h = 0)
dev.off()


