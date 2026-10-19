# Simulator for the space-time ETAS model with the Gaussian spatial kernel (mver = 2) of the
# ETAS package, conditioned on the observed doublet and on the observed events outside the
# target window. Coordinates are the projected coordinates stored in fit$object$revents.
# Copyright (C) 2026  Davood Kazemi Lafmejani and Elham Shabani. GPL-3.0-or-later.
sim_etas_gauss <- function(fit, beta, W_free, fixed_mm = 2.5, max_events = 5000) {
  obj <- fit$object; rv <- obj$revents; th <- fit$param
  win <- obj$region.win; T1 <- obj$rtperiod[2]; T0 <- obj$rtperiod[1]
  inside <- function(x, y) spatstat.geom::inside.owin(x, y, win)
  xr <- range(win$bdry[[1]]$x); yr <- range(win$bdry[[1]]$y)
  # events kept as observed: events outside the window and the (two) mainshocks
  fixed <- rv[, "flag"] == 0 | rv[, "mm"] >= fixed_mm
  F <- rv[fixed, c("tt", "xx", "yy", "mm"), drop = FALSE]
  rtgr <- function(n) { u <- runif(n); -log(1 - u * (1 - exp(-beta * W_free))) / beta }
  # background: Poisson number (expected integ0), uniform in time, spatial density = mixture of
  # Gaussian kernels centred at the observed events with weights = background probability
  nb <- rpois(1, th["mu"] * fit$integ0); B <- NULL   # expected background count is mu * integ0 (checked against ETAS source: mu * bk[j], mu * integ0)
  if (nb > 0) {
    got <- 0; xs <- ys <- numeric(0)
    while (got < nb) {
      i <- sample.int(nrow(rv), 1, prob = rv[, "prob"])
      x <- rv[i, "xx"] + rnorm(1, 0, fit$bwd[i]); y <- rv[i, "yy"] + rnorm(1, 0, fit$bwd[i])
      if (inside(x, y)) { xs <- c(xs, x); ys <- c(ys, y); got <- got + 1 }
    }
    B <- cbind(tt = runif(nb, T0, T1), xx = xs, yy = ys, mm = rtgr(nb))
  }
  all_free <- B; parents <- rbind(F, B)
  repeat {
    if (is.null(parents) || nrow(parents) == 0) break
    n <- rpois(nrow(parents), th["A"] * exp(th["alpha"] * parents[, "mm"]))
    if (sum(n) == 0) break
    i <- rep(seq_len(nrow(parents)), n)
    t <- parents[i, "tt"] + th["c"] * (runif(length(i))^(-1 / (th["p"] - 1)) - 1)
    s <- th["D"] * exp(th["gamma"] * parents[i, "mm"])   # D*exp(gamma*m) is the standard deviation (ETAS ffun2)
    x <- parents[i, "xx"] + rnorm(length(i)) * s; y <- parents[i, "yy"] + rnorm(length(i)) * s
    keep <- t <= T1 & inside(x, y)
    if (!any(keep)) break
    ch <- cbind(tt = t[keep], xx = x[keep], yy = y[keep], mm = rtgr(sum(keep)))
    all_free <- rbind(all_free, ch); parents <- ch
    if (!is.null(all_free) && nrow(all_free) > max_events) return(NULL)
  }
  ev <- rbind(F, all_free); ev <- ev[order(ev[, "tt"]), , drop = FALSE]
  flag <- as.integer(inside(ev[, "xx"], ev[, "yy"]))
  newobj <- obj
  newobj$revents <- cbind(ev, flag = flag, bkgd = 0, prob = 1, lambd = 0)
  newobj
}
