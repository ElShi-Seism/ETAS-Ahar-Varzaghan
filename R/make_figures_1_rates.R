suppressMessages(library(ETAS))
fit <- readRDS("out/reference_fit.rds")
ob <- fit$object; print(range(ob$region.poly$long)); print(range(ob$region.poly$lat))
res <- suppressMessages(resid.etas(fit)); tau <- res$tau; tt <- ob$revents[ob$revents[,"flag"]==1,"tt"]
cat(length(tau), length(tt), "\n"); print(head(tau)); print(tail(tau,2))
# rates over a map slightly larger than the window
lr <- c(38.20,38.65); gr <- c(46.5,47.0)
rt <- rates(fit, lat.range=lr, long.range=gr, plot.it=FALSE)
saveRDS(list(rt=rt, poly=ob$region.poly, tau=tau, tt=tt), "figdata.rds")
