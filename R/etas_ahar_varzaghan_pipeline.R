# ETAS Analysis Pipeline for the Ahar-Varzaghan 2012 Doublet Sequence
# Space-time ETAS fitting, residual analysis and changepoint detection,
# built on the R packages ETAS and changepoint.
# Copyright (C) 2026  Davood Kazemi Lafmejani and Elham Shabani
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ============================================================================
# ETAS analysis pipeline — Ahar-Varzaghan 2012 earthquake doublet sequence
# ============================================================================
# Dependencies: ETAS (Jalilian, 2019), changepoint (Killick & Eckley, 2014)
# License: GPL-3 (see LICENSE); depends on ETAS, which is GPL (>= 2).
# See README.md for background and CHANGELOG.md for a list of methodological
# fixes applied relative to an earlier version of this analysis.
# ============================================================================

library(ETAS)
library(changepoint)

SEED <- 20250923
ROUNDOFF <- TRUE            # ETAS::catalog() default: adds uniform jitter to coordinates
MIN_EVENTS_FOR_ETAS <- 50   # minimum events required to fit an independent ETAS model
STUDY_START <- as.Date("2012-08-11")
STUDY_END   <- as.Date("2024-11-26")
STUDY_DAYS  <- as.numeric(STUDY_END - STUDY_START)

# NOTE ON RANDOMNESS: ETAS::catalog() (roundoff = TRUE) adds uniform noise to
# epicentre coordinates BEFORE applying the lat/long window, so the number of
# target events near the window edge can change between unseeded runs.
# The seed is therefore set immediately before every catalog() call.
set.seed(SEED)

# ----------------------------------------------------------------------------
# Step 0: load and validate the raw catalog
# ----------------------------------------------------------------------------
catalog_path <- "ascii_correct.txt"
stopifnot(file.exists(catalog_path))

raw <- read.table(catalog_path, header = FALSE, sep = "\t",
                   col.names = c("long", "lat", "year", "month", "day",
                                 "mag", "sec", "hour", "minute"))

stopifnot(all(raw$hour >= 0 & raw$hour <= 23))
stopifnot(all(raw$minute >= 0 & raw$minute <= 59))
stopifnot(length(unique(raw$hour)) > 5)   # guards against a corrupted/constant hour field
stopifnot(!anyNA(raw))
cat("Raw events loaded:", nrow(raw), " | hour range:", range(raw$hour), "\n")

date_iso <- sprintf("%04d-%02d-%02d", raw$year, raw$month, raw$day)
time_iso <- sprintf("%02d:%02d:%05.2f", raw$hour, raw$minute, raw$sec)
k252 <- data.frame(date = date_iso, time = time_iso,
                    long = raw$long, lat = raw$lat, mag = raw$mag)

# ----------------------------------------------------------------------------
# Step 1: build the catalog object
# ----------------------------------------------------------------------------
kk252.cat <- catalog(k252,
                      time.begin  = as.character(STUDY_START),
                      study.start = as.character(STUDY_START),
                      study.end   = as.character(STUDY_END),
                      lat.range   = c(38.35, 38.5),
                      long.range  = c(46.5, 47),
                      mag.threshold = 3.2, roundoff = ROUNDOFF)
print(kk252.cat)

# ----------------------------------------------------------------------------
# Step 2: fit the base model with the Gaussian spatial kernel (mver = 2)
# ----------------------------------------------------------------------------
# param0 must have length 8 even under mver = 2 (the q slot is ignored but
# must be positive, per the package's internal validation).
param0_g <- c(mu = 0.5, A = 0.2, c = 0.05, alpha = 1.5,
              p = 1.2, D = 0.02, q = 1, gamma = 0.03)

fit0 <- etas(kk252.cat, param0 = param0_g, mver = 2)
cat("\n=== Fitted parameters (mver = 2) ===\n")
print(fit0$param)
ll_M0 <- max(fit0$loglik[fit0$loglik != 0])
cat("log-likelihood:", ll_M0, "\n")

# ----------------------------------------------------------------------------
# Step 2b: stability check with independent random starting points
# ----------------------------------------------------------------------------
stability_results <- data.frame()
for (i in 1:3) {
  rp <- c(mu = runif(1, 0.1, 2), A = runif(1, 0.05, 0.6), c = runif(1, 0.01, 0.15),
          alpha = runif(1, 0.8, 2.2), p = runif(1, 1.0, 1.5), D = runif(1, 0.001, 0.05),
          q = 1, gamma = runif(1, 0.1, 1.5))
  f <- tryCatch(etas(kk252.cat, param0 = rp, mver = 2), error = function(e) NULL)
  if (!is.null(f)) {
    ll <- max(f$loglik[f$loglik != 0])
    stability_results <- rbind(stability_results, data.frame(run = i, t(f$param), loglik = ll))
  }
}
cat("\n=== Stability check (3 independent random starts) ===\n")
print(stability_results)
d_spread <- diff(range(stability_results$D)) / mean(stability_results$D)
cat(sprintf("Relative spread of D across runs: %.4f\n", d_spread))
if (d_spread > 0.1) {
  warning("Parameters are still unstable even under mver = 2; further investigation needed.")
}

# ----------------------------------------------------------------------------
# Step 2c: comparison fit with the default power-law kernel (mver = 1)
# ----------------------------------------------------------------------------
# Included for transparency: this documents the identifiability problem that
# motivated the switch to mver = 2 above.
param0_p <- c(mu = 0.5, A = 0.2, c = 0.05, alpha = 1.5, p = 1.2, D = 0.02, q = 2.3, gamma = 0.03)
fit0_mver1 <- tryCatch(etas(kk252.cat, param0 = param0_p, mver = 1), error = function(e) NULL)
if (!is.null(fit0_mver1)) {
  cat("\n=== For comparison: mver = 1 (power-law kernel) result ===\n")
  print(fit0_mver1$param)
  cat(sprintf(">>> mver = 1 gives q = %.2f (a very large q means the power-law kernel is near its\n", fit0_mver1$param["q"]))
  cat(">>> Gaussian limit; D and q are then weakly identified and run-dependent).\n")
  cat(">>> mver = 2 is used as the primary reported model.\n")
}

# ----------------------------------------------------------------------------
# Step 3: residual analysis using the correct per-event transformed time
# ----------------------------------------------------------------------------
# NOTE: use res$tau, NOT res$tres. The latter is a binned diagnostic residual
# on an artificial fixed-length grid (n.temp - 1 points), not the classical
# per-event Ogata transformed time.
res <- resid.etas(fit0)
tau <- res$tau
cat("\ntau length:", length(tau), " | monotonically non-decreasing:", all(diff(tau) >= 0), "\n")
print(ks.test(res$U, punif))

# ----------------------------------------------------------------------------
# Step 4: structural breakpoint detection (PELT) with a sensitivity check
# ----------------------------------------------------------------------------
delta_tau <- diff(tau)
cpt_result <- cpt.mean(delta_tau, method = "PELT", penalty = "BIC", minseglen = 20)
cpt_idx <- cpts(cpt_result)
cat("\nBreakpoint(s) (event index, BIC penalty):", cpt_idx, "\n")

cat("Sensitivity to minseglen: ")
for (ms in c(10, 15, 20, 25)) {
  cp <- cpt.mean(delta_tau, method = "PELT", penalty = "BIC", minseglen = ms)
  cat(sprintf("[minseglen=%d -> %s] ", ms, paste(cpts(cp), collapse = ",")))
}
cat("\n")

# Map breakpoint index to a real calendar date using the catalog's own
# recorded elapsed times — NOT a proportional (uniform-spacing) estimate.
flg <- fit0$object$revents[, "flag"]
tt_target <- fit0$object$revents[flg == 1, "tt"]
if (length(cpt_idx) > 0) {
  for (idx in cpt_idx) {
    real_day <- tt_target[idx]
    real_date <- as.Date("2012-08-11") + real_day
    cat(sprintf("index %d -> day %.2f -> date %s\n", idx, real_day, as.character(real_date)))
  }
}

# ----------------------------------------------------------------------------
# Step 5: piecewise (before/after breakpoint) model comparison
# ----------------------------------------------------------------------------
# If a segment has fewer than MIN_EVENTS_FOR_ETAS events, an independent
# 8-parameter ETAS fit is not attempted; a simple observed rate is reported
# instead, and the formal likelihood-ratio comparison is skipped.
if (length(cpt_idx) > 0) {
  break_day <- tt_target[cpt_idx[1]]
  k252$elapsed <- as.numeric(difftime(
    as.POSIXct(paste(k252$date, k252$time), tz = "GMT"),
    as.POSIXct("2012-08-11", tz = "GMT"), units = "days"))

  in_win <- k252$lat >= 38.35 & k252$lat <= 38.5 &
            k252$long >= 46.5 & k252$long <= 46.99999999 & k252$mag >= 3.2
  seg1 <- k252[k252$elapsed <= break_day, c("date", "time", "long", "lat", "mag")]
  seg2 <- k252[k252$elapsed >  break_day, c("date", "time", "long", "lat", "mag")]
  n1 <- sum(in_win & k252$elapsed <= break_day)
  n2 <- sum(in_win & k252$elapsed >  break_day)
  cat(sprintf("\nSegment 1 target events: %d | Segment 2 target events: %d (window and Mc applied)\n", n1, n2))

  set.seed(SEED)
  cat_seg1 <- catalog(seg1, time.begin = as.character(STUDY_START),
                       study.start = as.character(STUDY_START),
                       study.end = as.character(STUDY_START + break_day),
                       lat.range = c(38.35, 38.5), long.range = c(46.5, 47),
                       mag.threshold = 3.2, roundoff = ROUNDOFF)
  fit_seg1 <- tryCatch(etas(cat_seg1, param0 = param0_g, mver = 2), error = function(e) NULL)

  fit_seg2 <- NULL
  if (n2 < MIN_EVENTS_FOR_ETAS) {
    seg2_days <- STUDY_DAYS - break_day
    cat(sprintf(">>> Segment 2 has only %d events (below the %d-event threshold);\n", n2, MIN_EVENTS_FOR_ETAS))
    cat(">>> no independent ETAS fit is attempted. Descriptive comparison instead:\n")
    cat(sprintf("Segment 2 observed rate: %.5f events/day (%d events over %.0f days)\n",
                n2 / seg2_days, n2, seg2_days))
    # Expected count after the break under the stationary model M0: the
    # compensator increment tau[n] - tau[idx] over the same interval. Observed
    # count in the open interval (t_idx, t_n) is n - idx - 1.
    # Post-selection caveat: the break was chosen from the same residuals, so
    # this p-value is descriptive, not a formal test.
    n_total <- length(tau)
    idx1 <- cpt_idx[1]
    expected2 <- tau[n_total] - tau[idx1]
    observed2 <- n_total - idx1 - 1
    p_pois <- ppois(observed2, expected2)
    cat(sprintf("Expected count after break under M0: %.1f | observed: %d | lower-tail Poisson p = %.3g (descriptive)\n",
                expected2, observed2, p_pois))
  } else {
    set.seed(SEED)
    cat_seg2 <- catalog(seg2, time.begin = as.character(STUDY_START + break_day),
                         study.start = as.character(STUDY_START + break_day),
                         study.end = as.character(STUDY_END),
                         lat.range = c(38.35, 38.5), long.range = c(46.5, 47),
                         mag.threshold = 3.2, roundoff = ROUNDOFF)
    fit_seg2 <- tryCatch(etas(cat_seg2, param0 = param0_g, mver = 2), error = function(e) NULL)
  }

  if (!is.null(fit_seg1) && !is.null(fit_seg2)) {
    ll1 <- max(fit_seg1$loglik[fit_seg1$loglik != 0])
    ll2 <- max(fit_seg2$loglik[fit_seg2$loglik != 0])
    ll_M1 <- ll1 + ll2
    df_LRT <- length(fit0$param) - 1   # q is not estimated under mver = 2
    LRT <- 2 * (ll_M1 - ll_M0)
    p_val <- pchisq(LRT, df = df_LRT, lower.tail = FALSE)
    cat(sprintf("\n=== M0 (stationary) vs M1 (piecewise) ===\nLRT=%.3f  df=%d  p=%.5f\n", LRT, df_LRT, p_val))
  } else {
    cat("\n>>> A formal M0/M1 likelihood-ratio comparison is not possible due to\n")
    cat(">>> insufficient data in segment 2.\n")
  }
}

# ----------------------------------------------------------------------------
# Save outputs
# ----------------------------------------------------------------------------
dir.create("ETAS_Results", showWarnings = FALSE)
saveRDS(list(fit0 = fit0, stability = stability_results, res = res, cpt = cpt_result,
             tt_target = tt_target, seed_used = SEED, n_target = length(tau), sessionInfo = sessionInfo()),
        file = "ETAS_Results/analysis_workspace.rds")

cat("\n=== Pipeline complete ===\n")
