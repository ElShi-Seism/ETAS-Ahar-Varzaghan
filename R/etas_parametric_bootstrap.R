# Parametric bootstrap for the Gaussian-kernel ETAS fit of the Ahar-Varzaghan sequence.
# For each replicate: (1) simulate a catalog from the fitted model (etas_sim_gauss.R), conditioned on
# the observed doublet and on the observed events outside the target window; (2) refit the model with
# ETAS::etas(mver = 2) starting at the fitted values; (3) recompute the diagnostics of the main pipeline
# on the refitted model. The null distributions therefore include the uncertainty of the estimated
# parameters, unlike the time-rescaling Monte Carlo of Step 4b.
#
# Usage (from the directory containing ascii_correct.txt):
#   Rscript R/etas_parametric_bootstrap.R run  <worker> <n_workers> <B> <outdir>   # simulate + refit
#   Rscript R/etas_parametric_bootstrap.R summarise <outdir>                        # tables from the CSVs
# Copyright (C) 2026  Davood Kazemi Lafmejani and Elham Shabani. GPL-3.0-or-later.
args <- commandArgs(trailingOnly = TRUE)
SEED <- 20250923
MC <- 3.2; CUT <- as.numeric(Sys.getenv("MMAX_FREE", "5.4")); W_FREE <- CUT - MC     # free events: magnitudes from the Gutenberg-Richter law truncated at 5.4,
                                  # the largest observed magnitude other than the doublet
here <- function(f) if (file.exists(f)) f else file.path("R", f)
source(here("etas_helpers.R")); source(here("etas_sim_gauss.R"))

obs_stats <- function(fit) {
  res <- resid.etas(fit); tau <- res$tau; dtau <- diff(tau)
  sc <- scan_exp(dtau, 20)
  ks <- suppressWarnings(ks.test(res$U, punif))
  c(ks_D = unname(ks$statistic), scan_two = unname(sc["two"]), split_two = unname(sc["k_two"]),
    scan_one = unname(sc["one"]), split_one = unname(sc["k_one"]), dev = max_dev_tau(tau),
    lag1 = stats::acf(dtau, lag.max = 1, plot = FALSE)$acf[2], lb10 = lb_stat(dtau, 10),
    runs_absz = abs(runs_z(dtau)),
    scan10 = unname(scan_exp(dtau, 10)["two"]), scan25 = unname(scan_exp(dtau, 25)["two"]))
}

if (args[1] == "run") {
  worker <- as.integer(args[2]); nw <- as.integer(args[3]); B <- as.integer(args[4]); outdir <- args[5]
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  library(ETAS)
  ref <- readRDS(file.path(outdir, "reference_fit.rds"))     # fit of the observed catalog (Mc = 3.2)
  beta <- unname(branching_ratio(ref, ref$object$revents[ref$object$revents[, "flag"] == 1, "mm"] + MC, MC)$beta)
  th <- ref$param; p0 <- c(th[1:6], q = 1, th[7]); names(p0) <- c("mu", "A", "c", "alpha", "p", "D", "q", "gamma")
  outfile <- file.path(outdir, sprintf("boot_worker%d.csv", worker))
  done <- if (file.exists(outfile)) read.csv(outfile)$b else integer(0)
  for (b in seq(worker, B, by = nw)) {
    if (b %in% done) next
    set.seed(SEED + 2000 + b)
    failfile <- file.path(outdir, sprintf("boot_failed_worker%d.csv", worker))
    note <- function(status, n_t = NA, n_tot = NA)
      write.table(data.frame(b = b, status = status, n_target = n_t, n_total = n_tot), failfile, sep = ",",
                  row.names = FALSE, col.names = !file.exists(failfile), append = file.exists(failfile))
    sim <- sim_etas_gauss(ref, beta, W_FREE)
    if (is.null(sim)) { note("simulation_exploded"); next }
    n_t <- sum(sim$revents[, "flag"] == 1); t0 <- Sys.time(); fit <- NULL; st <- NULL
    invisible(capture.output(fit <- tryCatch(suppressWarnings(etas(sim, param0 = p0, mver = 2)), error = function(e) NULL)))
    if (!is.null(fit)) invisible(capture.output(st <- tryCatch(obs_stats(fit), error = function(e) NULL)))
    if (is.null(st)) { note("fit_failed", n_t, nrow(sim$revents)); next }
    br <- tryCatch(branching_ratio(fit, sim$revents[sim$revents[, "flag"] == 1, "mm"] + MC, MC, Mmax = CUT), error = function(e) NULL)
    row <- data.frame(b = b, n_target = n_t, n_total = nrow(sim$revents),
                      loglik = max(fit$loglik[fit$loglik != 0]), t(fit$param),
                      t(setNames(fit$asd[fit$itr, ], paste0("se_", names(fit$param)))), t(st),
                      n_branch_cut = if (is.null(br)) NA else unname(br$truncated["estimate"]),
                      secs = as.numeric(difftime(Sys.time(), t0, units = "secs")), row.names = NULL)
    write.table(row, outfile, sep = ",", row.names = FALSE, col.names = !file.exists(outfile), append = file.exists(outfile))
  }
  quit(save = "no")
}

if (args[1] == "summarise") {
  outdir <- args[2]; library(ETAS)
  ref <- readRDS(file.path(outdir, "reference_fit.rds"))
  invisible(capture.output(obs <- obs_stats(ref)))
  fs <- list.files(outdir, "^boot_worker[0-9]+\\.csv$", full.names = TRUE)
  d <- do.call(rbind, lapply(fs, read.csv)); d <- d[order(d$b), ]
  ff <- list.files(outdir, "^boot_failed_worker[0-9]+\\.csv$", full.names = TRUE)
  nf <- if (length(ff)) nrow(do.call(rbind, lapply(ff, read.csv))) else 0
  B <- nrow(d); cat(sprintf("Replicates: %d completed, %d failed or exploded\n", B, nf))
  cat(sprintf("Simulated target events: median %.0f (5%%-95%%: %.0f-%.0f); observed 251\n", median(d$n_target),
              quantile(d$n_target, .05), quantile(d$n_target, .95)))
  pn <- c("mu", "A", "c", "alpha", "p", "D", "gamma"); th <- ref$param; se0 <- ref$asd[ref$itr, ]
  tab <- do.call(rbind, lapply(pn, function(k) {
    x <- d[[k]]; sx <- d[[paste0("se_", k)]]
    data.frame(parameter = k, estimate = th[k], asymptotic_se = se0[k], boot_mean = mean(x), boot_sd = sd(x),
               q025 = quantile(x, .025), q975 = quantile(x, .975),
               coverage_asymptotic = mean(abs(x - th[k]) <= 1.96 * sx), row.names = NULL)
  }))
  cat("\nParameter distribution over refitted simulated catalogs (coverage = share of replicates whose asymptotic 95% interval contains the generating value):\n")
  print(tab, digits = 4, row.names = FALSE)
  pv <- function(sim, o) (1 + sum(sim >= o)) / (length(sim) + 1)
  ptab <- data.frame(
    statistic = c("scan LR, two-sided (minseg 20)", "scan LR, one-sided (minseg 20)", "scan LR (minseg 10)", "scan LR (minseg 25)",
                  "max |N(tau)-tau| (early excess)", "K-S D", "Ljung-Box lag 10", "|lag-1 autocorrelation|", "runs test |z|"),
    observed = c(obs["scan_two"], obs["scan_one"], obs["scan10"], obs["scan25"], obs["dev"], obs["ks_D"], obs["lb10"], abs(obs["lag1"]), obs["runs_absz"]),
    p_bootstrap = c(pv(d$scan_two, obs["scan_two"]), pv(d$scan_one, obs["scan_one"]), pv(d$scan10, obs["scan10"]), pv(d$scan25, obs["scan25"]),
                    pv(d$dev, obs["dev"]), pv(d$ks_D, obs["ks_D"]), pv(d$lb10, obs["lb10"]), pv(abs(d$lag1), abs(obs["lag1"])), pv(d$runs_absz, obs["runs_absz"])),
    null_q95 = c(quantile(d$scan_two, .95), quantile(d$scan_one, .95), quantile(d$scan10, .95), quantile(d$scan25, .95), quantile(d$dev, .95),
                 quantile(d$ks_D, .95), quantile(d$lb10, .95), quantile(abs(d$lag1), .95), quantile(d$runs_absz, .95)))
  cat("\nBootstrap-calibrated p-values (the observed statistics are recomputed from the reference fit):\n"); print(ptab, digits = 4, row.names = FALSE)
  cat(sprintf("\nBranching ratio (cutoff 5.4) of the refits: mean %.2f, 95%% range %.2f-%.2f\n", mean(d$n_branch_cut, na.rm = TRUE),
              quantile(d$n_branch_cut, .025, na.rm = TRUE), quantile(d$n_branch_cut, .975, na.rm = TRUE)))
  write.csv(tab, file.path(outdir, "bootstrap_parameters.csv"), row.names = FALSE)
  write.csv(ptab, file.path(outdir, "bootstrap_pvalues.csv"), row.names = FALSE)
}
