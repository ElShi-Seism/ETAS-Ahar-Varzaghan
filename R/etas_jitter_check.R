# Sensitivity of the Gaussian-kernel fit to the coordinate jitter of ETAS::catalog (roundoff) and to the seed.
# Needs ascii_correct.txt in the working directory. Copyright (C) 2026 D. Kazemi Lafmejani and E. Shabani. GPL-3.0-or-later.
library(ETAS)
STUDY_START <- "2012-08-11"; STUDY_END <- "2024-11-26"
raw <- read.table("ascii_correct.txt", sep="\t", col.names=c("long","lat","year","month","day","mag","sec","hour","minute"))
k <- data.frame(date=sprintf("%04d-%02d-%02d",raw$year,raw$month,raw$day), time=sprintf("%02d:%02d:%05.2f",raw$hour,raw$minute,raw$sec), long=raw$long, lat=raw$lat, mag=raw$mag)
p0 <- c(mu=0.5,A=0.2,c=0.05,alpha=1.5,p=1.2,D=0.02,q=1,gamma=0.03)
cfgs <- list(list(r=FALSE,s=1), list(r=TRUE,s=101), list(r=TRUE,s=102), list(r=TRUE,s=103), list(r=TRUE,s=104))
out <- NULL
for (cf in cfgs) {
  set.seed(cf$s)
  cat_ <- catalog(k, time.begin=STUDY_START, study.start=STUDY_START, study.end=STUDY_END, lat.range=c(38.35,38.5), long.range=c(46.5,47), mag.threshold=3.2, roundoff=cf$r)
  f <- tryCatch(suppressWarnings(etas(cat_, param0=p0, mver=2)), error=function(e) NULL)
  if (is.null(f)) next
  out <- rbind(out, data.frame(roundoff=cf$r, seed=cf$s, n=sum(cat_$revents[,"flag"]==1), ll=max(f$loglik[f$loglik!=0]), t(f$param)))
  print(tail(out,1)); write.csv(out,"jitter.csv",row.names=FALSE)
}
