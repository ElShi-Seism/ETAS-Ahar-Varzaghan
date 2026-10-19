suppressMessages(library(ETAS)); source("etas_ownlik.R")
f <- readRDS("out/reference_fit.rds"); rv <- f$object$revents; th <- f$param
win <- f$object$region.win; xr <- range(win$bdry[[1]]$x); yr <- range(win$bdry[[1]]$y); T1 <- f$object$rtperiod[2]
d <- make_data(rv, f$integ0, T1, xr, yr)
nm <- names(th)
# observed information in natural parameters
Hn <- optimHess(th, function(p) -ll_etas(p, d), control = list(ndeps = th * 1e-4))
se_own <- sqrt(diag(solve(Hn))); cat("Hessian SE (own likelihood, natural scale):\n"); print(round(setNames(se_own, nm), 5))
# profile likelihood intervals
prof <- function(k, grid) {
  out <- sapply(grid, function(v) {
    free <- setdiff(seq_along(th), k)
    fn <- function(lp) { p <- th; p[free] <- exp(lp); p[k] <- v; -ll_etas(p, d) }
    o <- optim(log(th[free]), fn, method = "BFGS", control = list(maxit = 200, reltol = 1e-9)); -o$value })
  data.frame(value = grid, ll = out, drop = 2 * (f$loglikfv[f$itr] - out))
}
ci <- function(pr) { ok <- pr$drop <= 3.84; c(lo = min(pr$value[ok]), hi = max(pr$value[ok])) }
res <- list()
for (k in c(5, 6, 4, 3, 7, 2)) {
  v <- th[k]; grid <- v * c(seq(0.5, 0.9, by = 0.1), 0.95, 1, 1.05, seq(1.1, 2, by = 0.1))
  if (k %in% c(4, 5)) grid <- v + seq(-0.45, 0.45, by = 0.03) * v / 1.2 * 0.5
  if (k == 5) grid <- seq(1.05, 1.40, by = 0.02)
  pr <- prof(k, grid); res[[nm[k]]] <- pr; cat(nm[k], "profile 95% interval (LR=3.84):", signif(ci(pr), 4), " | grid", signif(range(grid), 3), "\n")
}
saveRDS(res, "out/profile_results.rds")
