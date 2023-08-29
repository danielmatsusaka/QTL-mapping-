version
install.packages("qtl")
install.packages("ASMap")


install.packages("gtools")
install.packages("fields")
install.packages("spam")
install.packages("dotCall64")
install.packages("grid")
install.packages("RColorBrewer")
install.packages("lattice")

library(qtl)
library(ASMap)
library(gtools)
library(fields)
library(spam)
library(dotCall64)
library(grid)
library(RColorBrewer)
library(lattice)

## Introduction

#The ColxPat genetic map Diego made is 10x longer than expected.  He used the MSTmap online tool.  I am going to try using the R/ASMap package, which uses the same algorithms, but allows for more flexibility and troubleshooting.

Will just try Chr1 first.

xxx <- pullCross(chr1.dat, type = "missing", pars = list(miss.thresh =
                                                            0.20))
xxx$missing$table


## 1. Load data

### Change input file to do other chromosomes

up.dat <- read.table("ch4_mst.txt",stringsAsFactors=FALSE)
#View(up.dat)
#coloque 2 en ves de 0 map size, cambie as.cross a false
chr1.dat <- mstmap.data.frame(up.dat, pop.type = "RIL6", dist.fun = "kosambi", objective.fun = "COUNT", p.value = 1e-06, noMap.dist = 5, noMap.size = 2, miss.thresh = 1, mvest.bc =FALSE, as.cross = TRUE, return.imputed = FALSE, trace = TRUE)

#up.dat <- read.table("C:/Users/Computer/Desktop/foto_paper/mapa en R/ch1/PARA_VERCH1.csv")


#2. Assess and remove missing data
plotMissing(chr1.dat)
sg <- statGen(chr1.dat, bychr = FALSE, stat.type = "miss")
hist(sg$miss)
#let's only keep lines with <10% missing genotypes
#cambie a 25% era 0.1
chr1.dat <- subset(chr1.dat, ind = sg$miss < (ncol(up.dat)*0.1)) #140 individuals
plotMissing(chr1.dat)

#missing data
#help("profileMark")
profileMark(chr1.dat, stat.type = c("seg.dist", "prop", "miss","dxo"), crit.val =
              "bonf", layout = c(1, 8), type = "l", cex = 0.5)
# remove markers with high proportion of missing data  # colocale 0.25
chr1.dat <- pullCross(chr1.dat, type = "missing", pars = list(miss.thresh =
                                                                0.25))

# remove markers with extreme segregation distortion
gt <- geno.table(chr1.dat)
gt[-log10(gt$P.value) > 25,] # this cutoff is arbitrary - there is a lot of segregation distortion here and i only want to remove this very badly behaving marker!
#--- puede que ponerle quiza ayude

todrop <- rownames(gt[-log10(gt$P.value) > 25,])

chr1.dat <- drop.markers(chr1.dat, todrop)
profileMark(chr1.dat, stat.type = c("seg.dist", "prop", "miss","dxo"), crit.val =
              "bonf", layout = c(1, 8), type = "l", cex = 0.5)
#3make a map for chromsome 1  # baja el p.value a -6
mapCHR1 <- mstmap(chr1.dat, bychr = FALSE, trace = TRUE, dist.fun ="kosambi", p.value = 1e-6)

chrlen(mapCHR1)

heatMap(mapCHR1)


## check number of crossovers per line
pg <- profileGen(mapCHR1, bychr = FALSE, stat.type = c("xo", "dxo","miss"), id = "Genotype", xo.lambda = 14, layout = c(1, 3), lty = 2, cex =0.7)

## two lines are excessive.
test <- subsetCross(mapCHR1, ind = !pg$xo.lambda)
mapTEST <- mstmap(test, bychr = TRUE, dist.fun = "kosambi", trace = TRUE,p.value = 1e-6) # bajale a -6
chrlen(mapTEST)

test.pg <- profileGen(mapTEST, bychr = FALSE, stat.type = c("xo", "dxo","miss"), id = "Genotype", xo.lambda = 14, layout = c(1, 3), lty = 2, cex =0.7)

hist(test.pg$stat$xo)

#4. Try re-estimating genetic distances

mapTESTe <- quickEst(mapTEST)
chrlen(mapTEST)
chrlen(mapTESTe)
#View(mapTESTe)

#write.csv(mapTESTe$geno,file="chromosoma1.csv")

ch1<-mapTESTe[["geno"]][["L"]][["map"]]
ch1

write.table(ch1,file = "ch5_MST_FINAL_miss_thresh_25.txt")
#write.table(mapTESTe, file = "cromo1.txt", sep = "\t",
          #  row.names = TRUE)
