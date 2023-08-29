
library(qtl)


df<-read.cross("csv",dir="C:/Users/Computer/Desktop",file="ctmr.csv",estimate.map=FALSE,na.strings = c("-",".","NA"),genotypes = c("A","B"),map.function = c("kosambi"),crosstype="riself")

summary(df)

help("summary.cross")
plotMissing(df)
par(mfrow=c(1,2), las=1)
plot(ntyped(df), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(df, "mar"), ylab="No. typed individuals", main="No. genotypes by marker")


help("ntyped")
plot(nmissing(df))
df <- subset(df, ind=(ntyped(df)>50))
help("subset")
nt.bymar <- ntyped(df, "mar")
nt.bymar
todrop <- names(nt.bymar[nt.bymar < 100])
df <- drop.markers(df, todrop)
totmar(df)

cg <- comparegeno(df)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])

wh<- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh



g <- pull.geno(df)
table(g[1,], g[15,])

#table(g[214,], g[216,])
#table(g[238,], g[288,])




!#for(i in 1:nrow(wh)) {tozero <- !is.na(g[wh[i,1],]) & !is.na(g[wh[i,2],]) & g[wh[i,1],] != g[wh[i,2],]df$geno[[1]]$data[wh[i,1],tozero] <- NA}
  
  #df <- subset(df, ind=-wh[,2])
  #print(dup <- findDupMarkers(df, exact.only=FALSE))
  
  #findDupMarkers(df)
  drop.dupmarkers(df)
print(dup <- findDupMarkers(df, exact.only=FALSE))
totmar(df)
#Look for markers with distorted segregation patterns
gt <- geno.table(df)
gt[gt$P.value < 0.05/totmar(df),]

View((gt))
#todrop <- rownames(gt[gt$P.value < 1e-10,])
todrop <- rownames(gt[gt$P.value < 1e-40,])
df <- drop.markers(df, todrop)

#Study individuals' genotype frequencies
g <- pull.geno(df)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)
for(i in 2:3)
  plot(gfreq[i,], ylab="Genotype frequency", main=c("A", "B")[i],ylim=c(0,1))





df<-markerlrt(df)
df<- est.rf(df)
checkAlleles(df, threshold=5)


rf <- pull.rf(df)
lod <- pull.rf(df, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

lg <- formLinkageGroups(df, max.rf=0.35, min.lod=6)
table(lg[,2])
nmar(df)
df <- formLinkageGroups(df, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)
plotRF(df, alternate.chrid=TRUE)


rf <- pull.rf(df)
lod <- pull.rf(df, what="lod")
#mn4 <- markernames(df, chr=4)
mn4 <- markernames(df, chr=1)
par(mfrow=c(2,1))
plot(rf, mn4[3], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
abline(h=0.5, lty=2)
plot(lod, mn4[3], bandcol="gray70", alternate.chrid=TRUE)
geno.crosstab(df, mn4[3], mn4[1])

mn5 <- markernames(df, chr=1)
geno.crosstab(df, mn4[3], mn5[1])

#toswitch <- markernames(df, chr=c(5, 7:11))
#df <- switchAlleles(df, toswitch)

df <- est.rf(df)
plotRF(df, alternate.chrid=TRUE)




rf <- pull.rf(df)
lod <- pull.rf(df, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")


lg <- formLinkageGroups(df, max.rf=0.35, min.lod=6)
table(lg[,2])
#Form linkage groups
df <- formLinkageGroups(df, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)
plotRF(df)

write.cross(df, file = "df.csv")
help("write.cross")
#cromosomo5 
# para este caso el 1
df <- orderMarkers(df, chr=1)
pull.map(df, chr=1)

rip1 <- ripple(df, chr=1, window=7)

rip1lik <- ripple(df, chr=1, window=4, method="likelihood",error.prob=0.005)



compareorder(df, chr=5, c(1:7,9,8), error.prob=0.01)

compareorder(df, chr=5, c(1:7,9,8), error.prob=0.001)

compareorder(df, chr=5, c(1:7,9,8), error.prob=0)

df <- switch.order(df, chr=5, c(1:7,9,8), error.prob=0.005)
pull.map(df, chr=5)


#Order markers on chromosome 4
# tenemos error no hay cambio de los marcadores 
df <- orderMarkers(df, chr=4)
pull.map(df, chr=4)

rip4 <- ripple(df, chr=4, window=7)
summary(rip4)

rip4lik <- ripple(df, chr=4, window=4, method="likelihood",error.prob=0.005)
summary(rip4lik)
df<- switch.order(df, chr=4, c(1:8,10,9), error.prob=0.005)
pull.map(df, chr=4)
### cromosoma 3
df <- orderMarkers(df, chr=3)
pull.map(df, chr=3)

rip3 <- ripple(df, chr=3, window=7)
summary(rip3)

rip3 <- ripple(df, chr=3, window=7)
summary(rip3)

rip3lik <- ripple(df, chr=3, window=4, method="likelihood",error.prob=0.005)
summary(rip3lik)


#cromosoma2

df <- orderMarkers(df, chr=2)
pull.map(df, chr=2)
rip2 <- ripple(df, chr=2, window=7)
summary(rip2)

rip2lik <- ripple(df, chr=2, window=4, method="likelihood",error.prob=0.005)
summary(rip2lik)


#df <- switch.order(df, chr=2, c(1,2,4,3,5:24), error.prob=0.005)
#pull.map(df, chr=2)


#Order markers on chromosome 1

df <- orderMarkers(df, chr=1)
pull.map(df, chr=1)

rip1 <- ripple(df, chr=1, window=7)
summary(rip1)

rip1lik <- ripple(df, chr=1, window=4, method="likelihood",error.prob=0.005)

summary(rip1lik)


#para terminar
summaryMap(df)

plotMap(df, show.marker.names=TRUE)
plotRF(df)

messedup <- switch.order(df, chr=1, c(1:11,23:33,12:22),error.prob=0.005)
plotRF(messedup, chr=1)

plotMap(messedup, show.marker.names=TRUE)
#Drop one marker at a time
dropone <- droponemarker(df, error.prob=0.005)

par(mfrow=c(2,1))
plot(dropone, lod=1, ylim=c(-100,0))
plot(dropone, lod=2, ylab="Change in chromosome length")
summary(dropone, lod.column=2)
badmar <- rownames(summary(dropone, lod.column=2))[1:3]
df <- drop.markers(df, badmar)

newmap <- est.map(df, error.prob=0.005)
mapthis <- replace.map(df, newmap)
summaryMap(df)

#Look for problem individuals
plot(countXO(df), ylab="Number of crossovers")
mapthis <- subset(df, ind=(countXO(df) < 50))
##para 5
summary(rip <- ripple(df, chr=5, window=7))
mapthis <- switch.order(df, chr=5, c(1:7,9,8), error.prob=0.005)
pull.map(df, chr=5)

newmap <- est.map(df, error.prob=0.005)
df <- replace.map(df, newmap)
summaryMap(df)


#Estimate genotyping error rate

loglik <- err <- c(0.001, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02)
for(i in seq(along=err)) {cat(i, "of", length(err), "\n")
  tempmap <- est.map(df, error.prob=err[i])
  loglik[i] <- sum(sapply(tempmap, attr, "loglik")) }
lod <- (loglik - max(loglik))/log(10)



plot(err, lod, xlab="Genotyping error rate", xlim=c(0,0.02), ylab=expression(paste(log[10], " likelihood")))

#Look for genotyping errors

df <- calc.errorlod(df, error.prob=0.005)
print(toperr <- top.errorlod(df, cutoff=6))

plotGeno(df, chr=1, ind=toperr$id[toperr$chr==1],cutoff=6, include.xo=FALSE)

df.clean <- df
for(i in 1:nrow(toperr)) {
  chr <- toperr$chr[i]
  id <- toperr$id[i]
  mar <- toperr$marker[i]
  df.clean$geno[[chr]]$data[df$pheno$id==id, mar] <- NA}

#Revisit segregation distortion
gt <- geno.table(df, scanone.output=TRUE)
par(mfrow=c(2,1))
plot(gt, ylab=expression(paste(-log[10], " P-value")))
plot(gt, lod=3:5, ylab="Genotype frequency")
abline(h=c(0.25, 0.5), lty=2, col="gray")
plotMap(df, show.marker.names=TRUE)


