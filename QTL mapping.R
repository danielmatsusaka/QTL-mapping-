library(qtl)
library(eqtl)
apa_Flowering<-read.cross("csv",dir="C:/Users/Computer/Desktop",file=".csv",na.strings = c("-","NA","."),genotypes = c("A","B"),alleles=c("A","B"),map.function = c("kosambi"),crosstype="riself")

nind(mapa_Flowering)
nphe(mapa_Flowering)
nchr(mapa_Flowering)
totmar(mapa_Flowering)
nmar(mapa_Flowering)



#Ajuste dos marcadores sobrepostos jittermap
Flowering_jm_cim <- jittermap(mapa_Flowering)
summary(Flowering_jm_cim)
plot.map(Flowering_jm_cim,show.marker.names=TRUE,horizontal=TRUE,shift=TRUE,alternate.chrid=FALSE)

plot.map(Flowering_jm_cim)
nmar(mapa_Flowering)
Flowering_jm_cim <- drop.nullmarkers(Flowering_jm_cim) 
totmar(Flowering_jm_cim)
Flowering_jm_cim <- est.rf(Flowering_jm_cim) 


Flowering_cim<- calc.genoprob(Flowering_jm_cim, step = 1, error.prob=0.0001)



#Flowering_cim Flowering_perm_cim


####flot1  no se encontraron qtls####

#N?vel de signific?ncia utilizando o m?todo de Haley-Knott

Flowering_perm_cim <- cim(Flowering_cim, method = "hk", pheno.col=1, n.perm=100000)
threshold_hk_cim <- summary(Flowering_perm_cim, alpha = 0.05)


plot(Flowering_perm_cim)
abline(v=threshold_hk_cim, lty = 2, col = "orange")


Flowering_qtl_cim = cim(Flowering_cim, method = "hk", pheno.col = 2,n.marcovar=3, map.function = "kosambi")


Flowering_perm_cim <- cim(Flowering_cim, method = "hk", pheno.col=2, n.perm=100000)
threshold_hk_cim <- summary(Flowering_perm_cim, alpha = 0.05)


plot(Flowering_perm_cim)
abline(v=threshold_hk_cim, lty = 2, col = "orange")


Flowering_qtl_cim = cim(Flowering_cim, method = "hk", pheno.col = 1,n.marcovar=3, map.function = "kosambi")

#View(Flowering_qtl_cim)
#write.csv(Flowering_qtl_cim, file = "Flowering_qtl_cimfv-1.csv")



sum_qtl_cim<-summary(Flowering_qtl_cim, perms = Flowering_perm_cim, alpha = 0.05, pvalues = TRUE)
#sum_qtl_cim<-summary(Flowering_qtl_cim, perms = Flowering_perm_cim, pvalues = TRUE)
sum_qtl_cim
plot(Flowering_qtl_cim, col = "blue")
plot(Flowering_qtl_cim, col = "blue",chr=c(1,2,3,4,5))
add.cim.covar(Flowering_qtl_cim, chr=c(1,2,3,4,5))
attr(Flowering_qtl_cim, "marker.covar.pos")
abline(h = threshold_hk_cim, lty = 2, col = "orange")
#Observando QTLs
qtl_cim <- makeqtl(Flowering_cim, chr=sum_qtl_cim$chr, pos=sum_qtl_cim$pos, what="prob")
plot(qtl_cim)

plot(Flowering_qtl_cim, chr=c(5), col =c("Blue"))

abline(h = threshold_hk_cim, lty = 2, col = "orange")

# grafico

myTicks<-c(0,1,2,3,4,5,6)
plot(Flowering_qtl_cim,lty=1,chr=c(1,2,3,4,5), col=c("blue"),ylab="LOD score", main="Phenotypes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylim=range(myTicks))
abline(h = threshold_hk_cim, lty = 2, col = "orange")

help("cex.axis")
# to find a interval mapping##
#The first and last rows define the ends of the intervals; the middle row is the estimated QTL location
lodint(Flowering_qtl_cim) 
bayesint(Flowering_qtl_cim) 
lodint(Flowering_qtl_cim, chr=5, expandtomarkers=TRUE)
bayesint(Flowering_qtl_cim, chr=5, expandtomarkers=TRUE)

lodint(Flowering_qtl_cim, chr=5, expandtomarkers=TRUE)
bayesint(Flowering_qtl_cim, chr=5, expandtomarkers=TRUE)

trait1 <- cim(Flowering_cim, pheno.col=1,
              window=15, method="hk", error.prob=0.001, 
              map.function=c("kosambi"))

trait1

#Observando os efeitos dos QTLs significativos nos cromossomos



Flowering_sim_cim <- sim.geno(Flowering_cim,n.draws=7, step=1, off.end=0, error.prob=0.0001,
                            map.function=c("kosambi"))
chr_all<-effectscan(Flowering_sim_cim, pheno.col=1, chr=c(1:5), get.se=TRUE, draw=TRUE,
                    gap=25, mtick=c("line","triangle"),add.legend=TRUE, alternate.chrid=FALSE, ylab="Effect", xlab="Linkage group")
out_cim <- fitqtl(Flowering_sim_cim, pheno.col=1, qtl=qtl_cim, method="hk", get.ests=TRUE)
summary(out_cim)
#para ch5 recuerda linkage chr[2] se?ala el qtl en orden ,igual position and lod
tablef1<-data.frame("Linkage Group"=c(sum_qtl_cim$chr[2]),
                    "Position"=c(sum_qtl_cim$pos[2]),
                    "LOD"=c(sum_qtl_cim$lod[2]),
                    "Additive effect"= c(summary(out_cim)$ests[2,1])
)
tablef1
#plot of genotype-specific phenotype means for 1 marker
mname1 <- find.marker(Flowering_cim, 5, 638.4432) # marker D1M437
effectplot(Flowering_cim, pheno.col=1, mname1=mname1)
### output of the function contains the means and SEs
output <- effectplot(Flowering_cim, mname1=mname1)
output

mname2 <- find.marker(Flowering_cim, 5,658.791 ) # marker D1M437
CH26.8<-effectplot(Flowering_cim, pheno.col=1, mname1=mname2)
CH26.8

output <- effectplot(Flowering_cim, mname1=mname2)
output
# PARA DOS MARCADORES
effectplot(Flowering_cim, mname1=mname1, mname2=mname2)

mname4 <- find.marker(Flowering_cim, 3, 181.022) # marker 94
CH94.0<-effectplot(Flowering_cim, pheno.col=5, mname1=mname4)
CH94.0

output1 <- effectplot(Flowering_cim, mname1=mname4)
output1

help("effectplot")

help("find.marker")


help("additive.effect")	
Flowering_cim<- calc.genoprob(Flowering_jm_cim, step = 1, error.prob=0.0001,off.end=0,map.function='kosambi', stepwidth='fixed')
help("sim.geno")

Flowering_sim_cim <- sim.geno( cross=Flowering_cim, step=1, error.prob=0.0001,off.end=0,map.function='kosambi', stepwidth='fixed')
## End(Not run)

# Genome scan and QTL detection
out.em <- scanone( Flowering_sim_cim, pheno.col=1:10, model='normal', method='hk')

View(out.em)
out.peak <- define.peak(out.em, 'all')
out.peak <- define.peak(out.em,lodcolumn='fitted_FLOWERING_LD',graph=TRUE,round=3)
#out.peak <- define.peak(out.em, 'all')
out.peak <- calc.adef(Flowering_sim_cim,out.em,out.peak,round=3)

out.peak
View(out.peak)

out.peak[[1]]$`5`$additive.effect
out.peak[[1]]$'5';


out.peak$fitted_FLOWERING_LD
out.peak$fitted_FLOWERING_LD$'5'
out.peak$fitted_FLOWERING_LD$'5'$additive.effect

#out.peak <- define.peak(out.em,lodcolumn=1,chr = 5,th=2.5,graph=FALSE,window.size=20,round=3, phe.name=1)
out.peak[["fitted_FLOWERING_LD"]]

out.peak[["fitted_WIGHT_LD"]]

#help(phe.names)
# Additive effect computing
out.peak <- calc.adef(Flowering_sim_cim1,out.em,out.peak)

out.peak






help("define.peak")
library(qtl)
#install.packages("eqtl")
library(eqtl)


#### flot2####

#N?vel de signific?ncia utilizando o m?todo de Haley-Knott


Flowering_perm_cim2 <- cim(Flowering_cim, method = "hk", pheno.col=2, n.perm=1000)
threshold_hk_cim <- summary(Flowering_perm_cim2, alpha = 0.05)

plot(Flowering_perm_cim2)
abline(v=threshold_hk_cim, lty = 2, col = "orange")



#Mapeamento CIM



Flowering_qtl_cim2 = cim(Flowering_cim, method = "hk", pheno.col = 2,n.marcovar=3, map.function = "kosambi")


sum_qtl_cim<-summary(Flowering_qtl_cim2, perms = Flowering_perm_cim2, alpha = 0.05, pvalues = TRUE)
sum_qtl_cim
plot(Flowering_qtl_cim2, col = "blue")
plot(Flowering_qtl_cim2, col = "blue",chr=c(1,2,3,4,5))
add.cim.covar(Flowering_qtl_cim2, chr=c(1,2,3,4,5))
attr(Flowering_qtl_cim2, "marker.covar.pos")
abline(h = threshold_hk_cim, lty = 2, col = "orange")


#Observando QTLs
qtl_cim <- makeqtl(Flowering_cim, chr=sum_qtl_cim$chr, pos=sum_qtl_cim$pos, what="prob")
plot(qtl_cim)

plot(Flowering_qtl_cim2, chr=c(5), show.marker.names=T, col =c("red"))
abline(h = threshold_hk_cim, lty = 2, col = "orange")

plot(Flowering_qtl_cim2, chr=c(5), col =c("lightblue"))


lodint(Flowering_qtl_cim2, chr=5, expandtomarkers=TRUE)
bayesint(Flowering_qtl_cim2, chr=5, expandtomarkers=TRUE)


#Observando os efeitos dos QTLs significativos nos cromossomos



chr_all<-effectscan(Flowering_sim_cim, pheno.col=2, chr=c(1:5), get.se=TRUE, draw=TRUE,
                    gap=25, mtick=c("line","triangle"),add.legend=TRUE, alternate.chrid=FALSE, ylab="Effect", xlab="Linkage group")
out_cim <- fitqtl(Flowering_sim_cim, pheno.col=2, qtl=qtl_cim, method="hk", get.ests=TRUE)

Flowering_sim_cim <- sim.geno(Flowering_cim,n.draws=7, step=1, off.end=0, error.prob=0.0001,
                            map.function=c("kosambi"))



summary(out_cim)
tablef2<-data.frame("Linkage Group"=c(sum_qtl_cim$chr[1]),
                    "Position"=c(sum_qtl_cim$pos[1]),
                    "LOD"=c(sum_qtl_cim$lod[1]),
                    "Additive effect"= c(summary(out_cim)$ests[2,1])
)
tablef2

tablef1<-data.frame("Linkage Group"=c(sum_qtl_cim$chr[2]),
                    "Position"=c(sum_qtl_cim$pos[2]),
                    "LOD"=c(sum_qtl_cim$lod[2]),
                    "Additive effect"= c(summary(out_cim)$ests[2,2])
)
tablef1

myTicks<-c(0,1,2,3,4,5,6)
plot(Flowering_qtl_cim2,lty=1,chr=c(1,2,3,4,5), col=c("blue"),ylab="LOD score", main="Phenotypes", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylim=range(myTicks))
abline(h = threshold_hk_cim, lty = 2, col = "orange")







##### para qtl luz


mname1 <- find.marker(Flowering_cim1, 2,83.405) # marker 260
effectplot(Flowering_cim1, pheno.col=1, mname1=mname1)
### output of the function contains the means and SEs
output <- effectplot(Flowering_cim1, mname1=mname1)
output



# PARA DOS MARCADORES



effectplot(Flowering_cim1, mname1=mname1, mname2=mname2)

mname4 <- find.marker(Flowering_cim1, 5, 638.443) # marker 260 y 267 iguales
CH26.0<-effectplot(Flowering_cim1, pheno.col=6, mname1=mname4)
CH26.0

output1 <- effectplot(Flowering_cim1, mname1=mname4)
output1

mname44 <- find.marker(Flowering_cim1, 3, 180.285) # marker 95 y 94 son iguales
CH9.4<-effectplot(Flowering_cim1, pheno.col=5, mname1=mname44)
CH9.4

mname5 <- find.marker(Flowering_cim1, 3, 160.025) # marker 79
CH9.4<-effectplot(Flowering_cim1, pheno.col=10, mname1=mname5)
CH9.4




help("effectplot")

help("find.marker")


help("additive.effect")	
Flowering_cim1<- calc.genoprob(Flowering_jm_cim1, step = 1, error.prob=0.0001,off.end=0,map.function='kosambi', stepwidth='fixed')
#help("sim.geno")

Flowering_sim_cim1 <- sim.geno( cross=Flowering_cim1, step=1, error.prob=0.0001,off.end=0,map.function='kosambi', stepwidth='fixed')
## End(Not run)

# Genome scan and QTL detection
out.em <- scanone( Flowering_sim_cim1, pheno.col=1:3, model='normal', method='hk')

View(out.em)
out.peak <- define.peak(out.em, 'all')
out.peak <- define.peak(out.em,lodcolumn='fitted_FLOWERING_LD',graph=TRUE,round=3)
#out.peak <- define.peak(out.em, 'all')
out.peak <- calc.adef(Flowering_sim_cim1,out.em,out.peak,round=3)

out.peak
View(out.peak)

out.peak[[1]]$`5`$additive.effect
out.peak[[1]]$'5';


out.peak$fitted_FLOWERING_LD
out.peak$fitted_FLOWERING_LD$'5'
out.peak$fitted_FLOWERING_LD$'5'$additive.effect

#out.peak <- define.peak(out.em,lodcolumn=1,chr = 5,th=2.5,graph=FALSE,window.size=20,round=3, phe.name=1)
out.peak[["fitted_FLOWERING_LD"]]

out.peak[["fitted_WIGHT_LD"]]

#help(phe.names)
# Additive effect computing
out.peak <- calc.adef(Flowering_sim_cim1,out.em,out.peak)

out.peak




