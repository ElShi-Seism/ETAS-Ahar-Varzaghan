# Sensitivity of the Ahar-Varzaghan ETAS results to the completeness magnitude Mc.
# Refits the Gaussian-kernel model (mver = 2) on the same window and period for
# several Mc values and repeats the diagnostics of the main pipeline. Run from the
# directory that contains ascii_correct.txt:  Rscript R/etas_sensitivity_Mc.R
# Copyright (C) 2026  Davood Kazemi Lafmejani and Elham Shabani. GPL-3.0-or-later.
library(ETAS)
library(changepoint)
SEED <- 20250923
MC_VALUES <- c(2.5, 2.7, 2.9, 3.2)
N_SIM <- 9999
HELPERS <- if (file.exists("etas_helpers.R")) "etas_helpers.R" else file.path("R", "etas_helpers.R")
stopifnot(file.exists(HELPERS)); source(HELPERS)

raw <- read.table("ascii_correct.txt", header = FALSE, sep = "\t",
                  col.names = c("long", "lat", "year", "month", "day", "mag", "sec", "hour", "minute"))
k252 <- data.frame(date = sprintf("%04d-%02d-%02d", raw$year, raw$month, raw$day),
                   time = sprintf("%02d:%02d:%05.2f", raw$hour, raw$minute, raw$sec),
                   long = raw$long, lat = raw$lat, mag = raw$mag)
param0_g <- c(mu = 0.5, A = 0.2, c = 0.05, alpha = 1.5, p = 1.2, D = 0.02, q = 1, gamma = 0.03)  # q: placeholder slot required by the package

rows <- list()
for (mc in MC_VALUES) {
  cat(sprintf("\n############ Mc = %.1f ############\n", mc))
  set.seed(SEED)   # same seed before every catalog() call
  ct <- catalog(k252, time.begin = "2012-08-11", study.start = "2012-08-11", study.end = "2024-11-26",
                lat.range = c(38.35, 38.5), long.range = c(46.5, 47), mag.threshold = mc, roundoff = TRUE)
  t0 <- Sys.time()
  fit <- etas(ct, param0 = param0_g, mver = 2)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  saveRDS(fit, sprintf("fit_Mc%.1f.rds", mc))
  rv <- fit$object$revents; flg <- rv[, "flag"] == 1; n <- sum(flg)
  ll <- max(fit$loglik[fit$loglik != 0]); se <- fit$asd[fit$itr, ]
  res <- resid.etas(fit); tau <- res$tau; dtau <- diff(tau)
  ks <- ks.test(res$U, punif)
  cp <- sapply(c(10, 15, 20, 25), function(ms) paste(cpts(cpt.mean(dtau, method = "PELT", penalty = "BIC", minseglen = ms)), collapse = "/"))
  tt <- rv[flg, "tt"]
  mcres <- mc_calibrate(tau, n_sim = N_SIM, minseg_set = c(10, 15, 20, 25), seed = SEED)
  ind <- residual_independence(dtau, n_sim = N_SIM, seed = SEED)
  br <- branching_ratio(fit, rv[flg, "mm"] + mc, Mc = mc)
  s20 <- mcres$scan[mcres$scan$minseglen == 20, ]
  k <- s20$split_two
  cat(sprintf("n = %d, loglik = %.2f, time %.0f s\n", n, ll, secs)); print(fit$param, digits = 4); print(se, digits = 3)
  cat("PELT breaks (minseglen 10/15/20/25):", cp, "\n"); print(mcres$scan, digits = 4)
  cat(sprintf("early excess: obs %.2f, q95 %.2f, p %.3f\n", mcres$dev_obs, mcres$dev_q95, mcres$dev_p))
  print(ind$lb, digits = 4); cat(sprintf("lag1 %.3f p %.4f | runs p %.4f\n", ind$lag1, ind$p_lag1, ind$p_runs))
  cat("branching (no cutoff):"); print(round(br$untruncated, 3)); cat("branching (cutoff):"); print(round(br$truncated, 3))
  rows[[length(rows) + 1]] <- data.frame(
    Mc = mc, n = n, loglik = ll, b = br$b,
    mu = fit$param["mu"], A = fit$param["A"], c = fit$param["c"], alpha = fit$param["alpha"],
    p = fit$param["p"], D = fit$param["D"], gamma = fit$param["gamma"],
    se_D = se["D"], se_alpha = se["alpha"], se_p = se["p"], se_c = se["c"],
    ks_p = ks$p.value, pelt20 = cp[3], split_scan20 = k, day_break = if (k > 0) tt[k] else NA,
    scan_stat = s20$stat_two, p_scan_two = s20$p_two, p_scan_one = s20$p_one,
    p_early = mcres$dev_p, p_lb10_mc = ind$lb$p_mc[ind$lb$lag == 10], p_runs = ind$p_runs,
    n_branch = br$untruncated["estimate"], n_branch_se = br$untruncated["se"],
    n_branch_cut = br$truncated["estimate"], n_branch_cut_se = br$truncated["se"], secs = secs,
    row.names = NULL)
  write.csv(do.call(rbind, rows), "sensitivity_Mc.csv", row.names = FALSE)
}
cat("\n=== Summary ===\n"); print(do.call(rbind, rows), digits = 4)
