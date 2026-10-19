# Independent implementation of the (conditional) ETAS log-likelihood of the ETAS package, mver = 2 (Gaussian
# kernel), for a rectangular window. Conditional on the background density u(x,y) and its integral integ0 stored by
# the package at its last declustering step:
#   ll = sum_{j target} log( mu*bk_j + sum_{i<j} kappa_i g(t_j-t_i) f(r_ji; m_i) )
#        - mu*integ0 - sum_{all j} kappa_j G(T-t_j) S_j
# with kappa = A exp(alpha m), g(t) = (p-1)/c (1+t/c)^-p, G(t) = 1-(1+t/c)^-(p-1), f Gaussian with standard
# deviation sigma = D exp(gamma m), S_j = integral of f over the window. Optional step change of the whole
# intensity by a factor r after time tb (model M1).
# Copyright (C) 2026  Davood Kazemi Lafmejani and Elham Shabani. GPL-3.0-or-later.
make_data <- function(rv, integ0, T1, xr, yr, bk = rv[, "bkgd"]) {
  N <- nrow(rv); t <- rv[, "tt"]; x <- rv[, "xx"]; y <- rv[, "yy"]
  dt <- outer(t, t, "-"); d2 <- outer(x, x, "-")^2 + outer(y, y, "-")^2
  list(N = N, t = t, x = x, y = y, m = rv[, "mm"], flag = rv[, "flag"] == 1, bk = bk, integ0 = integ0, T1 = T1,
       xr = xr, yr = yr, dt = dt, d2 = d2, low = dt > 0)
}
Gint <- function(s, c, p) ifelse(s > 0, 1 - (1 + pmax(s, 0) / c)^(-(p - 1)), 0)
ll_etas <- function(par, d, r = 1, tb = Inf) {
  mu <- par[1]; A <- par[2]; c <- par[3]; al <- par[4]; p <- par[5]; D <- par[6]; ga <- par[7]
  kap <- A * exp(al * d$m); sg <- D * exp(ga * d$m)
  g <- ifelse(d$low, (p - 1) / c * (1 + pmax(d$dt, 0) / c)^(-p), 0)
  S2 <- matrix(sg, d$N, d$N, byrow = TRUE)^2            # column i = parent
  f <- exp(-d$d2 / (2 * S2)) / (2 * pi * S2)
  trig <- rowSums(g * f * matrix(kap, d$N, d$N, byrow = TRUE))
  lam <- mu * d$bk + trig
  rj <- ifelse(d$t > tb, r, 1)
  s_j <- (pnorm((d$xr[2] - d$x) / sg) - pnorm((d$xr[1] - d$x) / sg)) * (pnorm((d$yr[2] - d$y) / sg) - pnorm((d$yr[1] - d$y) / sg))
  if (is.finite(tb)) {
    Gt <- Gint(tb - d$t, c, p); Gl <- Gint(d$T1 - d$t, c, p)
    comp_trig <- sum(kap * s_j * (Gt + r * (Gl - Gt)))
    comp_bg <- mu * d$integ0 * (tb / d$T1 + r * (1 - tb / d$T1))
  } else { comp_trig <- sum(kap * s_j * Gint(d$T1 - d$t, c, p)); comp_bg <- mu * d$integ0 }
  sum(log(pmax(rj * lam, 1e-25))[d$flag]) - comp_bg - comp_trig
}
fit_etas <- function(d, start, r_free_tb = NULL, r0 = 1, maxit = 300) {
  # parameters are optimised on the log scale; if r_free_tb is given, log r is an extra parameter
  if (is.null(r_free_tb)) {
    fn <- function(lp) -ll_etas(exp(lp), d)
    o <- optim(log(start), fn, method = "BFGS", control = list(maxit = maxit, reltol = 1e-10))
    list(par = exp(o$par), ll = -o$value, conv = o$convergence)
  } else {
    fn <- function(lp) -ll_etas(exp(lp[1:7]), d, r = exp(lp[8]), tb = r_free_tb)
    o <- optim(c(log(start), log(r0)), fn, method = "BFGS", control = list(maxit = maxit, reltol = 1e-10))
    list(par = exp(o$par[1:7]), r = exp(o$par[8]), ll = -o$value, conv = o$convergence)
  }
}
