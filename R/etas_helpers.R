# Helper functions for the Ahar-Varzaghan ETAS pipeline (GPL-3.0-or-later, see LICENSE).
# Sourced by etas_ahar_varzaghan_pipeline.R and etas_sensitivity_Mc.R.

# ---- MC BEGIN
scan_exp <- function(x, minseg = 20) {
  n <- length(x); S <- sum(x); cs <- cumsum(x)
  k <- minseg:(n - minseg)
  m1 <- cs[k] / k; m2 <- (S - cs[k]) / (n - k)
  llr <- 2 * (n * log(S / n) - k * log(m1) - (n - k) * log(m2))
  llr1 <- ifelse(m2 > m1, llr, 0)
  c(two = max(llr), k_two = k[which.max(llr)], one = max(llr1), k_one = k[which.max(llr1)])
}
max_dev_tau <- function(tau) {
  i <- seq_along(tau) - 1
  max(abs(i - (tau - tau[1])))
}
mc_calibrate <- function(tau, n_sim = 9999, minseg_set = c(10, 15, 20, 25), seed = 1) {
  set.seed(seed)
  dtau <- diff(tau); n <- length(dtau)
  obs_scan <- lapply(minseg_set, function(m) scan_exp(dtau, m))
  obs_dev <- max_dev_tau(tau)
  sim_two <- sim_one <- matrix(NA_real_, n_sim, length(minseg_set))
  sim_dev <- numeric(n_sim)
  for (b in seq_len(n_sim)) {
    x <- rexp(n)
    for (j in seq_along(minseg_set)) {
      s <- scan_exp(x, minseg_set[j]); sim_two[b, j] <- s["two"]; sim_one[b, j] <- s["one"]
    }
    sim_dev[b] <- max_dev_tau(c(0, cumsum(x)))
  }
  pv <- function(sim, obs) (1 + sum(sim >= obs)) / (length(sim) + 1)
  scan_tab <- data.frame(
    minseglen = minseg_set,
    split_two = sapply(obs_scan, function(s) s["k_two"]),
    stat_two  = sapply(obs_scan, function(s) s["two"]),
    p_two     = sapply(seq_along(minseg_set), function(j) pv(sim_two[, j], obs_scan[[j]]["two"])),
    split_one = sapply(obs_scan, function(s) s["k_one"]),
    stat_one  = sapply(obs_scan, function(s) s["one"]),
    p_one     = sapply(seq_along(minseg_set), function(j) pv(sim_one[, j], obs_scan[[j]]["one"])),
    row.names = NULL)
  list(scan = scan_tab, dev_obs = obs_dev, dev_p = pv(sim_dev, obs_dev),
       dev_q95 = unname(quantile(sim_dev, 0.95)), n_sim = n_sim)
}
# ---- MC END

# ---- Residual independence tests --------------------------------------------
# Under a correct model the transformed gaps are i.i.d. Exp(1). The K-S test ignores
# their order, so serial dependence is tested separately. Ljung-Box and the runs test
# are calibrated by Monte Carlo with i.i.d. Exp(1) data of the same length (the
# asymptotic chi-square law is not exact for exponential data). As for the scan test,
# the uncertainty of the estimated ETAS parameters is not propagated.
lb_stat <- function(x, lag) as.numeric(stats::Box.test(x, lag = lag, type = "Ljung-Box")$statistic)
runs_z <- function(x) {
  s <- x > stats::median(x); n1 <- sum(s); n2 <- sum(!s); r <- 1 + sum(s[-1] != s[-length(s)])
  mu <- 2 * n1 * n2 / (n1 + n2) + 1
  v <- 2 * n1 * n2 * (2 * n1 * n2 - n1 - n2) / ((n1 + n2)^2 * (n1 + n2 - 1))
  (r - mu) / sqrt(v)
}
residual_independence <- function(dtau, lags = c(5, 10, 20), n_sim = 9999, seed = 1) {
  set.seed(seed); n <- length(dtau)
  obs_lb <- sapply(lags, function(l) lb_stat(dtau, l)); obs_run <- abs(runs_z(dtau))
  obs_r1 <- stats::acf(dtau, lag.max = 1, plot = FALSE)$acf[2]
  sim_lb <- matrix(NA_real_, n_sim, length(lags)); sim_run <- numeric(n_sim); sim_r1 <- numeric(n_sim)
  for (b in seq_len(n_sim)) {
    x <- rexp(n)
    sim_lb[b, ] <- sapply(lags, function(l) lb_stat(x, l))
    sim_run[b] <- abs(runs_z(x)); sim_r1[b] <- stats::acf(x, lag.max = 1, plot = FALSE)$acf[2]
  }
  pv <- function(sim, obs) (1 + sum(sim >= obs)) / (length(sim) + 1)
  list(lag1 = obs_r1, p_lag1 = pv(abs(sim_r1), abs(obs_r1)),
       lb = data.frame(lag = lags, stat = obs_lb,
                       p_asymptotic = stats::pchisq(obs_lb, lags, lower.tail = FALSE),
                       p_mc = sapply(seq_along(lags), function(j) pv(sim_lb[, j], obs_lb[j]))),
       runs_absz = obs_run, p_runs = pv(sim_run, obs_run))
}

# ---- Branching ratio with delta-method standard error -----------------------
# n = A * beta / (beta - alpha) for magnitudes ~ exponential(beta) above Mc (no upper
# cutoff). With a cutoff at Mmax (width W = Mmax - Mc) the productivity is integrated
# over the truncated Gutenberg-Richter density. Var(beta) = beta^2 / N (Aki-Utsu), and
# beta is treated as independent of the ETAS parameters. Valid for mver = 2 only
# (parameter order mu, A, c, alpha, p, D, gamma; the covariance is corrected from fit$opt$avcov, see below).
branching_ratio <- function(fit, mags, Mc, Mmax = max(mags), dm = 0.1) {
  th <- fit$param
  # The package's fit$opt$avcov is inverted in scale (see README); the delta-method covariance of the
  # estimates is 16 * theta theta' * avcov, with theta the estimates (parameters = phi^2, phi = optimiser scale).
  V <- 16 * outer(th, th) * fit$opt$avcov
  N <- length(mags)
  b <- log10(exp(1)) / (mean(mags) - (Mc - dm / 2)); beta <- b * log(10)
  W <- Mmax - Mc
  nfun <- function(A, alpha, beta, W = Inf) {
    if (is.infinite(W)) return(A * beta / (beta - alpha))
    A * beta * (1 - exp(-(beta - alpha) * W)) / ((beta - alpha) * (1 - exp(-beta * W)))
  }
  one <- function(W) {
    f <- function(x) nfun(x[1], x[2], x[3], W)
    x0 <- c(th["A"], th["alpha"], beta); g <- numeric(3)
    for (i in 1:3) { h <- 1e-6 * max(1, abs(x0[i])); xp <- x0; xm <- x0; xp[i] <- xp[i] + h; xm[i] <- xm[i] - h
                     g[i] <- (f(xp) - f(xm)) / (2 * h) }
    S <- matrix(0, 3, 3)
    S[1:2, 1:2] <- V[c(2, 4), c(2, 4)]; S[3, 3] <- beta^2 / N
    est <- unname(f(x0)); se <- sqrt(drop(t(g) %*% S %*% g))
    c(estimate = est, se = se, lo = est - 1.96 * se, hi = est + 1.96 * se)
  }
  list(b = b, beta = beta, Mmax = Mmax, untruncated = one(Inf), truncated = one(W))
}
