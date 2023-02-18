#' Function to calculate Sig_eps, the covariance of a subject's latent z trajectory around the global group mean
#' @param len Squared exponential GP kernel length scale
#' @param sig_GP Squared exponential GP kernel variance
#' @param sig_RI Subject-specific random intercept variance
#' @param sig_resid Residual variance
Sig_eps_fn <- function(len, sig_GP, sig_RI, sig_resid) {
  Ones <- matrix(1, n_age, n_age)
  Iden <- diag(rep(1, n_age))
  Sig_GP <- outer(age_all, age_all, FUN = function(x, y) sig_GP * exp(- (x - y)^2 / 2 / len^2) )
  Sig_eps <- Sig_GP + Ones * sig_RI + Iden * sig_resid
  return(Sig_eps)
}
