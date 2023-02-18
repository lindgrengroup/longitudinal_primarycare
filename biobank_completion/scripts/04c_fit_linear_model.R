source("scripts/00_load_data.R")
list_import <- readRDS(file = paste0("output/linear_model_fitting/mat_vec_combine_task04a_", TASK_ID, ".RDS"))

Sig_beta <- Sig_eps_fn(len = 50, sig_GP = 1, sig_RI = 3, sig_resid = 0)
mu_beta <- rep(0, N_AGE_BIN)

n_param <- 1 + 2 * N_LOCI
beta_post_var <- solve(kronecker(X = diag(rep(1, n_param)), Y = solve(Sig_beta)) + list_import$mat)
beta_post_mn <- beta_post_var %*% list_import$vec

if (control$do_plots) {
  beta_gen <- mvtnorm::rmvnorm(n = 10, mean = mu_beta, sigma = Sig_beta)
  pdf("plots/test_beta_prior.pdf")
  matplot(t(beta_gen), ty = "l")
  dev.off()
}



