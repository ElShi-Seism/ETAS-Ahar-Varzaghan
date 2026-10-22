suppressMessages(library(ETAS))
fit <- readRDS("out/reference_fit.rds"); ob <- fit$object; x <- readRDS("figdata.rds")
rt <- x$rt
write.csv(rt$bkgd,"bkgd.csv",row.names=FALSE); write.csv(rt$clust,"clust.csv",row.names=FALSE)
write.csv(data.frame(x=rt$x),"gx.csv",row.names=FALSE); write.csv(data.frame(y=rt$y),"gy.csv",row.names=FALSE)
ev <- data.frame(long=ob$longlat$long, lat=ob$longlat$lat, mag=ob$revents[,"mm"]+3.2, flag=ob$revents[,"flag"], tt=ob$revents[,"tt"])
write.csv(ev,"events.csv",row.names=FALSE)
write.csv(data.frame(tt=x$tt,tau=x$tau),"tau.csv",row.names=FALSE)
print(table(ev$flag)); print(ev[ev$mag>=6,])
