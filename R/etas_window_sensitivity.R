# Sensitivity of the break to the extent of the spatial window (Gaussian kernel, Mc = 3.2).
# Usage (directory containing ascii_correct.txt and etas_helpers.R, i.e. run with R/ on the path):
#   Rscript R/etas_window_sensitivity.R [core_g mid_g full_g]
# Output: window_results.csv (see results/window_sensitivity.csv). Monte Carlo p: 2999 replicates, fixed parameters.
# Copyright (C) 2026  Davood Kazemi Lafmejani and Elham Shabani. GPL-3.0-or-later.
library(ETAS); source(if (file.exists("etas_helpers.R")) "etas_helpers.R" else "R/etas_helpers.R"); library(changepoint)
raw <- read.table("ascii_correct.txt", sep="\t", col.names=c("long","lat","year","month","day","mag","sec","hour","minute"))
cat_df <- data.frame(date=sprintf("%04d-%02d-%02d",raw$year,raw$month,raw$day), time=sprintf("%02d:%02d:%05.2f",raw$hour,raw$minute,raw$sec), long=raw$long, lat=raw$lat, mag=raw$mag)
# the slot named q is only a placeholder required by the package; the Gaussian kernel has no such parameter
p0 <- c(mu=0.5,A=0.2,c=0.05,alpha=1.5,p=1.2,D=0.02,q=1,gamma=0.03)
cfg <- list(
 core_g   = list(lat=c(38.35,38.5), lon=c(46.5,47.0), mver=2),
 mid_g    = list(lat=c(38.25,38.6), lon=c(46.4,47.0), mver=2),
 full_g   = list(lat=c(38.0,39.0),  lon=c(46.0,47.0), mver=2))
which <- commandArgs(TRUE); if (length(which)) cfg <- cfg[which]
out <- "window_results.csv"
for (nm in names(cfg)) {
  g <- cfg[[nm]]; t0 <- Sys.time(); set.seed(20250923)
  r <- tryCatch({
    ct <- catalog(cat_df, time.begin="2012-08-11", study.start="2012-08-11", study.end="2024-11-26", lat.range=g$lat, long.range=g$lon, mag.threshold=3.2, roundoff=TRUE)
    n <- sum(ct$revents[,"flag"]==1)
    invisible(capture.output(fit <- etas(ct, param0=p0, mver=g$mver)))
    ll <- max(fit$loglik[fit$loglik!=0]); k <- 7
    res <- resid.etas(fit); tau <- res$tau; dt <- diff(tau)
    sc <- scan_exp(dt, 20); mc <- mc_calibrate(tau, n_sim=2999, minseg_set=20, seed=1)
    tt <- ct$revents[ct$revents[,"flag"]==1,"tt"]
    kk <- as.integer(sc["k_two"]); N <- length(tau)
    pelt <- tryCatch(cpts(cpt.mean(dt, method="PELT", penalty="BIC", minseglen=20)), error=function(e) NA)
    ks <- suppressWarnings(ks.test(res$U, punif))
    exp_after <- tau[N]-tau[kk]; obs_after <- N-kk
    data.frame(config=nm, n_target=n, loglik=ll, AIC=-2*ll+2*k, t(fit$param), scan=unname(sc["two"]), split=kk, day=tt[kk+1], p_mc=mc$scan$p_two[1], pelt=paste(pelt,collapse=";"), obs_after=obs_after, exp_after=exp_after, ks_D=unname(ks$statistic), dev=max_dev_tau(tau), mean_before=mean(dt[1:kk]), mean_after=mean(dt[(kk+1):length(dt)]), mins=as.numeric(difftime(Sys.time(),t0,units="mins")))
  }, error=function(e) data.frame(config=nm, err=conditionMessage(e)))
  write.table(r, out, sep=",", row.names=FALSE, col.names=!file.exists(out), append=file.exists(out))
  print(r)
}
