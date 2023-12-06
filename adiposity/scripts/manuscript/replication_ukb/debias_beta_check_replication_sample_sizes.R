# Calculate debiasing factor

# For rs429358
b_obs <- 0.0597
se_obs <- 0.00485
afreq_obs <- 0.16
alpha_discovery <- 5E-8

z_stat <- b_obs/se_obs
z_crit <- qnorm(alpha_discovery, lower.tail = T)

debias_parameter <- se_obs*((dnorm(z_stat - z_crit) - dnorm(-z_stat - z_crit)) / 
                              (pnorm(z_stat - z_crit) + pnorm(-z_stat - z_crit)))

b_debiased <- b_obs - debias_parameter

# Estimate sample size needed for replication of debiased effect

# calculate non-centrality paramter for a range of N
getNCP <- function (sample_size, beta, afreq) {
  sigma <- sqrt(1 - 2*afreq*(1-afreq)*beta^2) 
  ses <- sigma/sqrt(sample_size*2*afreq*(1-afreq)) 
  return ((beta/ses)^2)
}

# calculate power for a given alpha and NCP
getPower <- function (alpha, ncp) {
  qthreshold <- qchisq(alpha, df = 1, ncp = 0, lower.tail = F)
  pwr <- pchisq(qthreshold, df = 1, ncp = ncp, lower = F)
  return (pwr)
}

alpha_rep_gws <- 5E-6
alpha_nominal <- 0.05

# Test a range of sample sizes

nrange <- seq(1, 1E5, by = 1E2)
ncp_for_nrange <- getNCP(nrange, beta = b_debiased, afreq = afreq_obs)
power_gws <- getPower(alpha = alpha_rep_gws,
                      ncp = ncp_for_nrange)

plot(nrange, power_gws)

# Sample size for 80% power
nrange[which.min(abs(power_gws - 0.8))]

# Calculate for all the discovered SNPs ----

VARIDS <- c("rs429358", "chr6_26076446", "rs11778922", "rs9467663", "rs61955499", "rs12953815")
b_obs <- c(0.0597,  0.011837299999999999, 0.0156858, 0.0107597, 0.06685039999999999, 0.0144239)
se_obs <- c(0.00485, 0.00197532, 0.00275548, 0.00178387, 0.0121087, 0.0025591999999999998)
afreq_obs <- c(0.156, 0.43678, 0.378749, 0.4185939999999999, 0.012036000000000047, 0.49634799999999996)

names(b_obs) <- VARIDS
names(se_obs) <- VARIDS
names(afreq_obs) <- VARIDS

alpha_discovery <- 5E-8
z_crit <- qnorm(alpha_discovery, lower.tail = T)

b_debiased <- lapply(VARIDS, function (vid) {
  z_stat <- b_obs[vid]/se_obs[vid]
  debias_parameter <- se_obs[vid]*((dnorm(z_stat - z_crit) - dnorm(-z_stat - z_crit)) / 
                                (pnorm(z_stat - z_crit) + pnorm(-z_stat - z_crit)))
  
  res <- b_obs[vid] + debias_parameter
  return (res)
})
b_debiased <- unlist(b_debiased)
names(b_debiased) <- VARIDS

# Estimate sample size needed for replication of debiased effect

# calculate non-centrality paramter for a range of N
getNCP <- function (sample_size, beta, afreq) {
  sigma <- sqrt(1 - 2*afreq*(1-afreq)*beta^2) 
  ses <- sigma/sqrt(sample_size*2*afreq*(1-afreq)) 
  return ((beta/ses)^2)
}

# calculate power for a given alpha and NCP
getPower <- function (alpha, ncp) {
  qthreshold <- qchisq(alpha, df = 1, ncp = 0, lower.tail = F)
  pwr <- pchisq(qthreshold, df = 1, ncp = ncp, lower = F)
  return (pwr)
}

alpha_rep_gws <- 0.05/6
alpha_nominal <- 0.05

# Test a range of sample sizes

nrange <- seq(1, 500000, by = 1000)

lapply(VARIDS, function (vid) {
  ncp_for_nrange <- getNCP(nrange, beta = b_debiased[vid], afreq = afreq_obs[vid])
  power_gws <- getPower(alpha = alpha_rep_gws,
                        ncp = ncp_for_nrange)
  
  plot(nrange, power_gws,
       main = paste0("rsid: ", vid))
  
  # Sample size for 80% power
  paste0("rsid: ", vid, " - ", nrange[which.min(abs(power_gws - 0.8))])
})


