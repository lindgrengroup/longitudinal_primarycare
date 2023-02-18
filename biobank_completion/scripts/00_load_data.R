##########################################
# R packages 
##########################################
require(argparse)
require(splines)
require(lubridate)
require(tidyverse)
require(RColorBrewer)
theme_set(theme_bw())

##################################################################
# Load functions, control parameters and data
##################################################################
for (fn_file in list.files("scripts/functions", full.names = TRUE)) {
  source(fn_file)
}
control <- get_control_parameters()

force_recreate_model_dat <- FALSE

##############################################################################
# Load data and split into list with each participant as an element
##############################################################################
dat <- readRDS(file.path(control$adiposity_root, "highdim_splines/data/dat_to_model.rds"))[[control$PHENO]][[control$SEX_STRATA]]
dat$age_days <- round(dat$age_t1 * 365.25 + dat$t_diff)
day_unique <- min(dat$age_days):max(dat$age_days)
N_DAY <- length(day_unique)
all_subj <- unique(dat$eid)
n_subj <- length(all_subj)
if (!"model_dat" %in% ls() | force_recreate_model_dat) {
  model_dat_file <- "output/model_dat.RDS"
  if (!file.exists(model_dat_file) | force_recreate_model_dat) {
    model_dat <- split(dat, f = dat$eid)
    saveRDS(model_dat, file = model_dat_file)
  } else {
    model_dat <- readRDS(file = model_dat_file)
  }
}
day_unique <- min(dat$age_days):max(dat$age_days)
N_DAY <- length(day_unique)
all_subj <- unique(dat$eid)
n_subj <- length(all_subj)


##############################################################################
# Load data and split into list with each participant as an element
##############################################################################
age_bin_df <- get_age_bin_df(control)
N_AGE_BIN <- nrow(age_bin_df)
subj_time_ranges <- data.frame(eid = all_subj)
dat <- dat[order(dat$age_days), ]
subj_time_ranges$lower <- dat[match(subj_time_ranges$eid, dat$eid), ] %>% pull("age_days")
dat <- dat[order(dat$age_days, decreasing = TRUE), ]
subj_time_ranges$upper <- dat[match(subj_time_ranges$eid, dat$eid), ] %>% pull("age_days")
age_bin_df$n_subj <- sapply(1:N_AGE_BIN, function(i) sum(subj_time_ranges$lower <= age_bin_df$upp_day[max(c(1, i - 1))] & 
                                                           subj_time_ranges$upper >= age_bin_df$low_day[min(c(N_AGE_BIN, i))]))

y_bin_df <- get_y_bin_df(control, y_vec = dat$value)
N_Y_BIN <- nrow(y_bin_df)
age_all <- age_bin_df$l
n_age <- length(age_all)

RANDOM_SEED <- 160522
set.seed(RANDOM_SEED)

##############################################################################
# Note current heterogeneous bin width
##############################################################################
if (control$do_plots) {
  pdf(file = "plots/bin_wid.pdf", width = 9, height = 12)
  plot(y_bin_df$u - y_bin_df$l, ty = "p", log = "y")
  dev.off()
}
# bin_centre_quantiles <- seq(1 / N_Y_BIN, 1, by = 1 / N_Y_BIN) - 1 / 2 / N_Y_BIN
# bin_centre_in <- quantile(dat$value, probs = bin_centre_quantiles)

##############################################################################
# Generate synthetic genetic data
##############################################################################
N_LOCI <- 100
snp_names <- paste0("SNP_", 1:N_LOCI)
genetic_data <- matrix(0, n_subj, N_LOCI, dimnames = list(all_subj, snp_names))
set.seed(1)
maf_vec <- runif(n = N_LOCI, min = 0, max = .5)
for (j in 1:N_LOCI) {
  genetic_data[, j] <- rbinom(n = n_subj, size = 2, prob = maf_vec[j])
}
sex_data <- as.numeric(runif(n = n_subj, min = 0, max = 1) < .5)

meta <- data.frame(subj_id = all_subj, sex = sex_data)
meta <- cbind(meta, genetic_data)
