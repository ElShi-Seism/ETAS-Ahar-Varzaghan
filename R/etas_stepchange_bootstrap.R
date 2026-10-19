# Competing model for "relative quiescence": the Gaussian-kernel ETAS intensity multiplied by a factor r after a
# change time tb (M1, 8 parameters) against the stationary model (M0, 7 parameters). tb is chosen by scanning
# (with the ETAS parameters fixed at the M0 estimate, closed form), then M1 is refitted with all parameters free.
# The likelihood-ratio statistic is calibrated by a parametric bootstrap that repeats the whole procedure
# (including the choice of tb) on catalogs simulated from M0 (etas_sim_gauss.R). The background density u(x,y) is
# held fixed at the reference estimate (the same in the data and in the simulations).
# Usage: Rscript etas_stepchange_bootstrap.R obs <outdir>
#        Rscript etas_stepchange_bootstrap.R run <worker> <n_workers> <B> <outdir> <seed_offset>
#        Rscript etas_stepchange_bootstrap.R summarise <outdir>
# Copyright (C) 2026  Davood Kazemi Lafmejani and Elham Shabani. GPL-3.0-or-later.
args <- commandArgs(trailingOnly = TRUE)
SEED <- 20250923; MC <- 3.2; W_FREE <- 5.4 - MC
here <- function(f) if (file.exists(f)) f else file.path("R", f)
suppressMessages({library(ETAS); library(Rcpp)})
source(here("etas_helpers.R")); source(here("etas_sim_gauss.R")); source(here("etas_ownlik.R")); source(here("etas_ownlik_cpp.R"))

setup <- function(outdir) {
  ref <- readRDS(file.path(outdir, "reference_fit.rds")); rv <- ref$object$revents
  win <- ref$object$region.win
  list(ref = ref, rv = rv, xr = range(win$bdry[[1]]$x), yr = range(win$bdry[[1]]$y), T1 = ref$object$rtperiod[2])
}
bk_ref <- function(S, x, y) {   # background density at (x, y): (1/T) sum_i pb_i dGauss(r_i, bwd_i), as in the package
  rv <- S$rv; w <- rv[, "prob"]; h <- S$ref$bwd
  sapply(seq_along(x), function(j) { r2 <- (x[j] - rv[, "xx"])^2 + (y[j] - rv[, "yy"])^2
    sum(w * exp(-r2 / (2 * h^2)) / (2 * pi * h^2)) / S$T1 })
}
lc <- function(par, D, r = 1, tb = Inf) ll_cpp(par, D$t, D$x, D$y, D$m, D$flag, D$bk, D$integ0, D$T1, D$xr, D$yr, r, tb)
to_u <- function(par) { u <- log(par); u[5] <- log(par[5] - 1); u }      # log-scale, with p > 1 enforced
from_u <- function(u) { par <- exp(u); par[5] <- 1 + exp(u[5]); par }
fitc <- function(D, start, tb = NULL, r0 = 1) {
  pen <- function(u) if (any(abs(u) > 12)) 1e6 else 0
  if (is.null(tb)) { fn <- function(u) { pp <- from_u(u); if (any(abs(u) > 12)) return(1e6); -lc(pp, D) }
    o <- optim(to_u(start), fn, method = "BFGS", control = list(maxit = 300, reltol = 1e-11))
    return(list(par = from_u(o$par), ll = -o$value, conv = o$convergence)) }
  fn <- function(u) { if (any(abs(u) > 12)) return(1e6); -lc(from_u(u[1:7]), D, exp(u[8]), tb) }
  o <- optim(c(to_u(start), log(r0)), fn, method = "BFGS", control = list(maxit = 300, reltol = 1e-11))
  list(par = from_u(o$par[1:7]), r = exp(o$par[8]), ll = -o$value, conv = o$convergence)
}
scan_tb <- function(par, D, minev = 20) {
  tt <- D$t[D$flag]; n <- length(tt); k <- minev:(n - minev); tbs <- tt[k]      # change right after the k-th target event
  kap <- par[2] * exp(par[4] * D$m); sg <- par[6] * exp(par[7] * D$m); c <- par[3]; p <- par[5]
  s_j <- (pnorm((D$xr[2] - D$x) / sg) - pnorm((D$xr[1] - D$x) / sg)) * (pnorm((D$yr[2] - D$y) / sg) - pnorm((D$yr[1] - D$y) / sg))
  Gl <- 1 - (1 + (D$T1 - D$t) / c)^(-(p - 1))
  Clate <- sapply(tbs, function(tb) { Gt <- ifelse(tb > D$t, 1 - (1 + pmax(tb - D$t, 0) / c)^(-(p - 1)), 0)
    par[1] * D$integ0 * (1 - tb / D$T1) + sum(kap * s_j * (Gl - Gt)) })
  nl <- n - k
  lr <- 2 * (nl * log(nl / Clate) - nl + Clate)
  list(tb = tbs[which.max(lr)], lr_cond = max(lr), r_cond = (nl / Clate)[which.max(lr)], k = k[which.max(lr)],
       lr_cond_def = max(ifelse(nl / Clate < 1, lr, 0)))
}
analyse <- function(D, start) {
  f0 <- fitc(D, start); sc <- scan_tb(f0$par, D)
  f1 <- fitc(D, f0$par, tb = sc$tb, r0 = sc$r_cond)
  c(ll0 = f0$ll, ll1 = f1$ll, LR_full = 2 * (f1$ll - f0$ll), LR_cond = sc$lr_cond, tb = sc$tb, k = sc$k, r1 = f1$r,
    r_cond = sc$r_cond, p0 = f0$par[5], p1 = f1$par[5], A0 = f0$par[2], A1 = f1$par[2], D0 = f0$par[6], conv = f0$conv + f1$conv)
}
if (args[1] == "obs") {
  S <- setup(args[2]); rv <- S$rv; th <- S$ref$par; if (is.null(th)) th <- S$ref$param
  D <- make_data(rv, S$ref$integ0, S$T1, S$xr, S$yr)
  D$bk <- rv[, "bkgd"]
  a <- analyse(D, th); print(a)
  cat(sprintf("AIC M0 %.2f, AIC M1 %.2f (conditional log-likelihoods)\n", -2 * a["ll0"] + 14, -2 * a["ll1"] + 16))
  saveRDS(a, file.path(args[2], "stepchange_observed.rds"))
}
if (args[1] == "run") {
  worker <- as.integer(args[2]); nw <- as.integer(args[3]); B <- as.integer(args[4]); outdir <- args[5]; off <- as.integer(args[6])
  S <- setup(outdir); th <- S$ref$param
  beta <- 0.860 * log(10)
  outfile <- file.path(outdir, sprintf("step_worker%d.csv", worker))
  done <- if (file.exists(outfile)) read.csv(outfile)$b else integer(0)
  for (b in seq(worker, B, by = nw)) {
    if (b %in% done) next
    set.seed(SEED + off + b)
    sim <- sim_etas_gauss(S$ref, beta, W_FREE); if (is.null(sim)) next
    rvs <- sim$revents; D <- make_data(rvs, S$ref$integ0, S$T1, S$xr, S$yr, bk = bk_ref(S, rvs[, "xx"], rvs[, "yy"]))
    a <- tryCatch(analyse(D, th), error = function(e) NULL); if (is.null(a)) next
    row <- data.frame(b = b, n_target = sum(rvs[, "flag"] == 1), t(a), row.names = NULL)
    write.table(row, outfile, sep = ",", row.names = FALSE, col.names = !file.exists(outfile), append = file.exists(outfile))
  }
}
if (args[1] == "summarise") {
  outdir <- args[2]; a <- readRDS(file.path(outdir, "stepchange_observed.rds"))
  d <- do.call(rbind, lapply(list.files(outdir, "^step_worker[0-9]+\\.csv$", full.names = TRUE), read.csv))
  pv <- function(s, o) (1 + sum(s >= o)) / (1 + length(s))
  cat(sprintf("Replicates: %d\nObserved: LR_full = %.2f, LR_cond = %.2f, r = %.3f at tb = %.1f (event %d)\n", nrow(d), a["LR_full"], a["LR_cond"], a["r1"], a["tb"], a["k"]))
  cat(sprintf("Bootstrap null of LR_full: median %.2f, 95%% quantile %.2f, 99%% %.2f; p = %.4f\n", median(d$LR_full), quantile(d$LR_full, .95), quantile(d$LR_full, .99), pv(d$LR_full, a["LR_full"])))
  cat(sprintf("Bootstrap null of LR_cond: 95%% quantile %.2f; p = %.4f\n", quantile(d$LR_cond, .95), pv(d$LR_cond, a["LR_cond"])))
  m <- d$n_target >= 200 & d$n_target <= 320
  cat(sprintf("Replicates with 200-320 target events: %d, p(LR_full) = %.4f\n", sum(m), pv(d$LR_full[m], a["LR_full"])))
  cat(sprintf("AIC difference observed: %.2f (M1 better if negative); share of null replicates with a lower AIC for M1: %.3f\n",
              2 * (8 - 7) - 2 * (a["ll1"] - a["ll0"]), mean(2 - 2 * (d$ll1 - d$ll0) < 0)))
}
