df<-read.csv(file.choose(),header = TRUE, na.strings = "NA",",")
df
View(df)
summary(df)
### INSTALL  PACKAGES ###
source("https://bioconductor.org/biocLite.R")
install.packages("Rcpp")
install.packages("stringi")
install.packages("reshape2")
install.packages("colorspace")
install.packages("ggplot2")
install.packages("corrplot")
install.packages("pbkrtest")
install.packages("quantreg")
install.packages("mgcv")
install.packages("devtools")
biocLite("VIM")
biocLite("zoo")
biocLite("sva")
biocLite("mvoutlier")
install.packages("sgeostat")
biocLite("ggbiplot")
install.packages("ggbiplot")
biocLite("lme4")
biocLite("MASS")



### LOAD LIBRARIES ###
require(XLConnect)
#install.packages("XLConnect")
#install.packages("XLConectJars")
library(XLConnect)
library(data.table)
library(reshape2)
library(ggplot2)
library(corrplot)
library(VIM)
library(zoo)
library(sva)
library(mvoutlier)
library(ggbiplot)
library(devtools)
#install.packages("processx")
library(lme4)
library(MASS)
library(lattice)
version
install.packages("lme4")
pheno=df
colnames(pheno)
class(pheno); dim(pheno); head(pheno)
### LOOK FOR MISSING DATA USING VIM LIBRARY
#aggr(pheno[c(7:1)])
#dat.missing = summary(aggr(pheno[c(7:7)], plot=TRUE))

### WE FIRST APPLY A BOXCOX TRANSFORMATION ###
### SEE HERE FOR MORE EXPLANATION : https://www.isixsigma.com/tools-templates/normality/making-data-normal-using-box-cox-power-transformation/

##### flowering TRIAT####

# run the box-cox transformation

bc.azul = boxcox(azul~1, data=pheno)

# We extract the correction lambda value
# The correction formulat depends on the lambda range value


## LdH: and it seems to be the maximum
lambda.azul = bc.azul$x[which.max(bc.azul$y)]


### LdH: See how normal is the data using a qqplot of the residuals
m.azul = lm(azul~1, data=pheno)

# Then we QQ-plot the residuals
#op <- par(pty = "s", mfrow = c(1, 2))
qqnorm(m.azul$residuals); qqline(m.azul$residuals)
# qqnorm(mnew.alphaT$residuals); qqline(mnew.alphaT$residuals)
#par(op)

# We use the shapiro test for checking normality of the residuals of the LM
# The null-hypothesis of this test is that the population is normally distributed. Thus, if the p-value is less than the chosen alpha level, then the null hypothesis is rejected and there is evidence that the data tested are not from a normally distributed population; in other words, the data are not normal. On the contrary, if the p-value is greater than the chosen alpha level, then the null hypothesis that the data came from a normally distributed population cannot be rejected (e.g., for an alpha level of 0.05, a data set with a p-value of 0.02 rejects the null hypothesis that the data are from a normally distributed population)
# LdH: If the pvalue > 0.05, the null hypothesis (normality) is not rejected
shapiro.m.azul = shapiro.test(m.azul$residuals)
par(mfrow=c(2,2))
hist(pheno$azul, main = "Raw data", xlab="azul", col ="grey60")
hist(sqrt(pheno$azul), main = "Sqrt Transformation", xlab="azul", col ="grey60")
hist(log10(pheno$azul), main = "Log10 Transformation", xlab="", col ="grey60")
hist(((pheno$azul)^lambda.azul), main = "BoxCox Transformation", xlab="azul", col ="grey60")

# We apply a final shapiro test on the data
shapiro.test(pheno$azul)
shapiro.test(sqrt(pheno$azul))
shapiro.test(log10(pheno$azul))
shapiro.test((pheno$azul)^lambda.azul)


# The BOXCOX Transformation shows the highest level of correction
# Thus, we transform the data according to 10^lambda
corr.data.azul=data.frame((pheno$azul)^lambda.azul)
colnames(corr.data.azul)=c("azul")
head(corr.data.azul)


## qqplots again to check the corrected data
m.azul = lm(azul~1, data=pheno)
m.corr.data.azul= lm(azul~1, data=corr.data.azul)
op <- par(pty = "s", mfrow = c(1, 2))
qqnorm(m.azul$residuals); qqline(m.azul$residuals)
qqnorm(m.corr.data.azul$residuals); qqline(m.corr.data.azul$residuals)

write.csv(corr.data.azul,file="corrdataluz_azul.csv")
colnames(pheno)

#### azul.osc TRAIT #### 

# run the box-cox transformation

bc.azul.osc = boxcox(azul.osc~1, data=pheno)
# We extract the correction lambda value
# The correction formulat depends on the lambda range value
lambda.azul.osc = bc.azul.osc$x[which.max(bc.azul.osc$y)]


m.azul.osc = lm(azul.osc~1, data=pheno)

# Then we QQ-plot the residuals
qqnorm(m.azul.osc$residuals); qqline(m.azul.osc$residuals)
# We use the shapiro test for checking normality of the residuals of the LM
shapiro.m.azul.osc = shapiro.test(m.azul.osc$residuals)

shapiro.m.azul.osc 

# WE PLOT THE COMPARISON BETWEEN DIFFERENT CORRECTIONS TO CROSSCHECK THE BOXCOX TRANSFORMATION
par(mfrow=c(2,2))
hist(pheno$azul.osc, main = "Raw data", xlab="azul.osc", col ="grey60")
hist(sqrt(pheno$azul.osc), main = "Sqrt Transformation", xlab="azul.osc", col ="grey60")
hist(log10(pheno$azul.osc), main = "Log10 Transformation", xlab="azul.osc", col ="grey60")
#hist(((pheno$azul.osc)^lambda.azul.osc), main = "BoxCox Transformation", xlab="azul.osc", col ="grey60")
hist(((pheno$azul.osc)^lambda.azul.osc), main = "BoxCox Transformation", xlab="azul.osc", col ="grey60")


# We apply a final shapiro test on the data
shapiro.test(pheno$azul.osc)
shapiro.test(sqrt(pheno$azul.osc))
shapiro.test(log10(pheno$azul.osc))
shapiro.test((pheno$azul.osc)^lambda.azul.osc)

# The BOXCOX Transformation shows the highest level of correction
# Thus, we transform the data according to 10^lambda
# We replace NA values by mean of the trait using the na.aggregate function
corr.data.azul.osc=data.frame((pheno$azul.osc)^lambda.azul.osc)
colnames(corr.data.azul.osc)=c("azul.osc")
head(corr.data.azul.osc)

### plot qqplot for before and after
m.azul.osc = lm(azul.osc~1, data=pheno)
m.corr.data.azul.osc = lm(azul.osc~1, data=corr.data.azul.osc)
op <- par(pty = "s", mfrow = c(1, 2))
qqnorm(m.azul.osc$residuals); qqline(m.azul.osc$residuals)
qqnorm(m.corr.data.azul.osc$residuals); qqline(m.corr.data.azul.osc$residuals)

write.csv(corr.data.azul.osc,file="corrdata_azul.osc.csv")

#### rojo TRAIT #### 
colnames(pheno)
# run the box-cox transformation
bc.rojo = boxcox(rojo~1, data=pheno)
# We extract the correction lambda value
# The correction formulat depends on the lambda range value

lambda.rojo = bc.rojo$x[which.max(bc.rojo$y)]

m.rojo = lm(rojo~1, data=pheno)

# Then we QQ-plot the residuals

qqnorm(m.rojo$residuals); qqline(m.rojo$residuals)


# We use the shapiro test for checking normality of the residuals of the LM
shapiro.m.rojo = shapiro.test(m.rojo$residuals)

shapiro.m.rojo

# WE PLOT THE COMPARISON BETWEEN DIFFERENT CORRECTIONS TO CROSSCHECK THE BOXCOX TRANSFORMATION
par(mfrow=c(2,2))
hist(pheno$rojo, main = "Raw data", xlab="rojo", col ="grey60")
hist(sqrt(pheno$rojo), main = "Sqrt Transformation", xlab="rojo", col ="grey60")
hist(log10(pheno$rojo), main = "Log10 Transformation", xlab="rojo", col ="grey60")
hist(((pheno$rojo)^lambda.rojo), main = "BoxCox Transformation", xlab="rojo", col ="grey60")

# We apply a final shapiro test on the data
shapiro.test(pheno$rojo)
shapiro.test(sqrt(pheno$rojo))
shapiro.test(log10(pheno$rojo))
shapiro.test((pheno$rojo)^lambda.rojo)

# The log10 Transformation shows the highest level of correction
# Thus, we transform the data according to log10

corr.data.rojo=data.frame((pheno$rojo)^lambda.rojo)
colnames(corr.data.rojo)=c("rojo")
head(corr.data.rojo)

### plot qqplot for before and after
m.rojo = lm(rojo~1, data=pheno)
m.corr.data.rojo = lm(rojo~1, data=corr.data.rojo)
op <- par(pty = "s", mfrow = c(1, 2))
qqnorm(m.rojo$residuals); qqline(m.rojo$residuals)
qqnorm(m.corr.data.rojo$residuals); qqline(m.corr.data.rojo$residuals)

write.csv(corr.data.rojo,file="corrdataluz_rojo.csv")

####rojo.osc TRAIT #### 
colnames(pheno)
# run the box-cox transformation
bc.rojo.osc = boxcox(rojo.osc~1, data=pheno)

# We extract the correction lambda value
# The correction formulat depends on the lambda range value

lambda.rojo.osc = bc.rojo.osc$x[which.max(bc.rojo.osc$y)]

m.rojo.osc = lm(rojo.osc~1, data=pheno)


# Then we QQ-plot the residuals
qqnorm(m.rojo.osc$residuals); qqline(m.rojo.osc$residuals)


# We use the shapiro test for checking normality of the residuals of the LM
shapiro.m.rojo.osc = shapiro.test(m.rojo.osc$residuals)

shapiro.m.rojo.osc

# WE PLOT THE COMPARISON BETWEEN DIFFERENT CORRECTIONS TO CROSSCHECK THE BOXCOX TRANSFORMATION
par(mfrow=c(2,2))
hist(pheno$rojo.osc, main = "Raw data", xlab="rojo.osc", col ="grey60")
hist(sqrt(pheno$rojo.osc), main = "Sqrt Transformation", xlab="rojo", col ="grey60")
hist(log10(pheno$rojo.osc), main = "Log10 Transformation", xlab="rojo", col ="grey60")
hist(((pheno$rojo.osc)^lambda.rojo.osc), main = "BoxCox Transformation", xlab="LENGHT", col ="grey60")

# We apply a final shapiro test on the data
shapiro.test(pheno$rojo.osc)
shapiro.test(sqrt(pheno$rojo.osc))
shapiro.test(log10(pheno$rojo.osc))
shapiro.test((pheno$rojo.osc)^lambda.rojo.osc)

# The BOXCOX Transformation shows the highest level of correction
corr.data.rojo.osc=data.frame((pheno$rojo.osc)^lambda.rojo.osc)
colnames(corr.data.rojo.osc)=c("rojo.osc")
head(corr.data.rojo.osc)




### plot qqplot for before and after
m.rojo.osc = lm(rojo.osc~1, data=pheno)
m.corr.data.rojo.osc = lm(rojo.osc~1, data=corr.data.rojo.osc)
op <- par(pty = "s", mfrow = c(1, 2))
qqnorm(m.rojo.osc$residuals); qqline(m.rojo.osc$residuals)
qqnorm(m.corr.data.rojo.osc$residuals);qqline(m.corr.data.rojo.osc$residuals)


write.csv(corr.data.rojo.osc,file="corrdataluz_rojo.osc.csv")
#### fr TRAIT #### 

colnames(pheno)
# run the box-cox transformation
bc.fr = boxcox(fr~1, data=pheno)
# We extract the correction lambda value
# The correction formulat depends on the lambda range value

lambda.fr = bc.fr$x[which.max(bc.fr$y)]

m.fr = lm(fr~1, data=pheno)

# Then we QQ-plot the residuals
qqnorm(m.fr$residuals); qqline(m.fr$residuals)

# We use the shapiro test for checking normality of the residuals of the LM
shapiro.m.fr = shapiro.test(m.fr$residuals)

shapiro.m.fr

# WE PLOT THE COMPARISON BETWEEN DIFFERENT CORRECTIONS TO CROSSCHECK THE BOXCOX TRANSFORMATION
par(mfrow=c(2,2))
hist(pheno$fr, main = "Raw data", xlab="deltaT", col ="grey60")
hist(sqrt(pheno$fr), main = "Sqrt Transformation", xlab="deltaT", col ="grey60")
hist(log10(pheno$fr), main = "Log10 Transformation", xlab="deltaT", col ="grey60")
hist(((pheno$fr)^lambda.fr), main = "BoxCox Transformation", xlab="deltaT", col ="grey60")

# We apply a final shapiro test on the data
shapiro.test(pheno$fr)
shapiro.test(sqrt(pheno$fr))
shapiro.test(log10(pheno$fr))
shapiro.test((pheno$fr)^lambda.fr)


# The log10 Transformation shows the highest level of correction
# Thus, we transform the data according to log10


corr.data.fr=data.frame((pheno$fr)^lambda.fr)
colnames(corr.data.fr)=c("fr")
head(corr.data.fr)




### plot qqplot for before and after




m.fr = lm(fr~1, data=pheno)
m.corr.data.fr = lm(fr~1, data=corr.data.fr)
op <- par(pty = "s", mfrow = c(1, 2))
qqnorm(m.fr$residuals); qqline(m.fr$residuals)
qqnorm(m.corr.data.fr$residuals);qqline(m.corr.data.fr$residuals)
write.csv(corr.data.fr,file="corrdataluz_fr.csv")





#### fr.osc TRAIT #### 

colnames(pheno)
# run the box-cox transformation
bc.fr.osc = boxcox(fr.osc~1, data=pheno)
# We extract the correction lambda value
# The correction formulat depends on the lambda range value

lambda.fr.osc = bc.fr.osc$x[which.max(bc.fr.osc$y)]

m.fr.osc = lm(fr.osc~1, data=pheno)

# Then we QQ-plot the residuals
qqnorm(m.fr.osc$residuals); qqline(m.fr.osc$residuals)

# We use the shapiro test for checking normality of the residuals of the LM
shapiro.m.fr.osc = shapiro.test(m.fr.osc$residuals)

shapiro.m.fr.osc

# WE PLOT THE COMPARISON BETWEEN DIFFERENT CORRECTIONS TO CROSSCHECK THE BOXCOX TRANSFORMATION
par(mfrow=c(2,2))
hist(pheno$fr.osc, main = "Raw data", xlab="deltaT", col ="grey60")
hist(sqrt(pheno$fr.osc), main = "Sqrt Transformation", xlab="deltaT", col ="grey60")
hist(log10(pheno$fr.osc), main = "Log10 Transformation", xlab="deltaT", col ="grey60")
hist(((pheno$fr.osc)^lambda.fr.osc), main = "BoxCox Transformation", xlab="fr.osc", col ="grey60")

# We apply a final shapiro test on the data
shapiro.test(pheno$fr.osc)
shapiro.test(sqrt(pheno$fr.osc))
shapiro.test(log10(pheno$fr.osc))
shapiro.test((pheno$fr.osc)^lambda.fr.osc)


# The log10 Transformation shows the highest level of correction
# Thus, we transform the data according to log10


corr.data.fr.osc=data.frame((pheno$fr.osc)^lambda.fr.osc)
colnames(corr.data.fr.osc)=c("fr.osc")
head(corr.data.fr.osc)




### plot qqplot for before and after




m.fr = lm(fr.osc~1, data=pheno)
m.corr.data.fr.osc = lm(fr.osc~1, data=corr.data.fr.osc)
op <- par(pty = "s", mfrow = c(1, 2))
qqnorm(m.fr.osc$residuals); qqline(m.fr.osc$residuals)
qqnorm(m.corr.data.fr.osc$residuals);qqline(m.corr.data.fr.osc$residuals)
write.csv(corr.data.fr.osc,file="corrdataluz_fr.osc.csv")




#### osc TRAIT #### 

colnames(pheno)
# run the box-cox transformation
bc.oscuridad = boxcox(oscuridad~1, data=pheno)
# We extract the correction lambda value
# The correction formulat depends on the lambda range value

lambda.oscuridad = bc.oscuridad$x[which.max(bc.oscuridad$y)]

m.oscuridad = lm(oscuridad~1, data=pheno)

# Then we QQ-plot the residuals
qqnorm(m.oscuridad$residuals); qqline(m.oscuridad$residuals)

# We use the shapiro test for checking normality of the residuals of the LM
shapiro.m.oscuridad = shapiro.test(m.oscuridad$residuals)

shapiro.m.oscuridad

# WE PLOT THE COMPARISON BETWEEN DIFFERENT CORRECTIONS TO CROSSCHECK THE BOXCOX TRANSFORMATION
par(mfrow=c(2,2))
hist(pheno$oscuridad, main = "Raw data", xlab="deltaT", col ="grey60")
hist(sqrt(pheno$oscuridad), main = "Sqrt Transformation", xlab="deltaT", col ="grey60")
hist(log10(pheno$oscuridad), main = "Log10 Transformation", xlab="deltaT", col ="grey60")
hist(((pheno$oscuridad)^lambda.oscuridad), main = "BoxCox Transformation", xlab="deltaT", col ="grey60")


# We apply a final shapiro test on the data
shapiro.test(pheno$oscuridad)
shapiro.test(sqrt(pheno$oscuridad))
shapiro.test(log10(pheno$oscuridad))
shapiro.test((pheno$oscuridad)^lambda.oscuridad)


# The log10 Transformation shows the highest level of correction
# Thus, we transform the data according to log10


corr.data.oscuridad=data.frame((pheno$oscuridad)^lambda.oscuridad)
colnames(corr.data.oscuridad)=c("oscuridad")
head(corr.data.oscuridad)




### plot qqplot for before and after




m.oscuridad = lm(oscuridad~1, data=pheno)
m.corr.data.oscuridad = lm(oscuridad~1, data=corr.data.oscuridad)
op <- par(pty = "s", mfrow = c(1, 2))
qqnorm(m.oscuridad$residuals); qqline(m.oscuridad$residuals)
qqnorm(m.corr.data.oscuridad$residuals);qqline(m.corr.data.oscuridad$residuals)
write.csv(corr.data.oscuridad,file="corrdataluz_oscuridad.csv")



### CREATE A DATAFRAME THAT CONTAINS THE TRANFORMED PHENOTYPIC DATA 
### BE CAREFUL WITH THE NUMBER OF COLUM TO EXTRACT FROM THE INITIAL 'pheno' DATAFRAME 

corr.df=cbind(pheno[,1:13], corr.data.alphaT, corr.data.aT3,
              corr.data.betaT, corr.data.gammaT,
              corr.data.deltaT, corr.data.Total_Toco)
dim(corr.df); head(corr.df)

corr.df=cbind(pheno[,1:7], corr.data.alphaT, corr.data.aT3,
              corr.data.betaT, corr.data.gammaT,
              corr.data.deltaT, corr.data.Total_Toco)
corr.df

df<-read.csv(file.choose(),header = TRUE,na.strings = "NA")
corr.df=
  # write.table(corr.df, "__corrected_phenotypic_dataframe.txt", sep="\t", col.names=TRUE, quote=FALSE)
  


### LOAD THE CORRECTED DATA
# load sd y ld data frame corr
df<-read.csv(file.choose(),header = TRUE,na.strings = "NA")
df=corr.df
corr.df
colnames(corr.df)
View(corr.df)
### USE WITH ALPHA TOCO DATA
### TEST FOR DIFFERENTS MODELS
### MORE TO READ ABOUT MODEL SELECTION : http://www.stat.umn.edu/geyer/5931/mle/sel.pdf
### about overfitting: We say it "overfits" the data, meaning it's too close to the data and not close enough to the true population regression function.

model0 = lmer(FLOWERING ~ (1|RIL), data=corr.df)
model0b = lmer(FLOWERING ~ RIL + (1|RIL), data = corr.df)

### add foto  (trataiento sd o ld)
model1 = lmer(FLOWERING ~ (1|RIL) + (1|foto), data=corr.df)
model1b = lmer(FLOWERING ~ RIL + (1|foto), data=corr.df)
model1c = lmer(FLOWERING ~ RIL + foto + (1|RIL:foto), data=corr.df)

### ADD EFFECT of 'mesured_in' parameter###
model1 = lmer(FLOWERING ~ (1|RIL) + (1|REP), data=corr.df)
model1b = lmer(FLOWERING ~ RIL + (1|REP), data=corr.df)
#model1c = lmer(FLOWERING ~ RIL + REP + (1|RIL:REP), data=corr.df)


# and its interaction with the genotype
#model2 = lmer(FLOWERING ~ (1|RIL) + (1|REP) + (1|RIL:REP), data=corr.df)
#model2b = lmer(FLOWERING ~ RIL + (1|RIL) + (1|REP) + (1|RIL:REP), data=corr.df)

### ADD THE EFFECT OF 'Treatment' parameter to model 1
#model3 = lmer(alphaT ~ (1|Genotype_name) + (1|Measured_in) + (1|Treatment), data=corr.df)
model3 = lmer(FLOWERING ~ (1|RIL) + (1|REP), data=corr.df)

# and its interaction with the genotype

#model4 = lmer(FLOWERING ~ (1|Genotype_name) + (1|Measured_in) + (1|Treatment) + (1|Genotype_name:Treatment), data=corr.df)

model4 = lmer(FLOWERING ~ (1|RIL) + (1|REP) + (1|foto), data=corr.df)

# add the interaction with the 'measured_in'
model4b = lmer(FLOWERING ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto), data=corr.df)
model5 = lmer(FLOWERING ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto) + (1|RIL:REP), data=corr.df)
#model6 = lmer(FLOWERING ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto) + (1|RIL:REP) + (1|foto:REP), data=corr.df)
#no le gusta el modelo6 modelo falla en coverngncia 

#model7 = lmer(FLOWERING ~ (1|RIL)  + (1|RIL:foto) + (1|RIL:REP), data=corr.df)
#modelo sobre ajustado


model8 = lmer(FLOWERING ~ (1|RIL) + (1|foto) + (1|RIL:foto) + (1|RIL:REP), data=corr.df)


anova(model0,model0b,model1,model1b,model3)


VarCorr(model0)
VarCorr(model0b)
VarCorr(model1)
VarCorr(model1b)
VarCorr(model3)

VarCorr(model1b)
VarCorr(model1c)
VarCorr(model2)
VarCorr(model2b)
VarCorr(model3)
VarCorr(model4)
VarCorr(model4b)
VarCorr(model5)
VarCorr(model8)


model1

ranef(model1)
rr1 = ranef(model1, condVar = TRUE)
plot1 = dotplot(rr1)[1]
plot2 = dotplot(rr1)[2]
plot3 = dotplot(rr1)[3]


print(VarCorr(model0),comp="Variance")
print(VarCorr(model1),comp="Variance")

fixef(model0)
fixef(model1)

ranef(model0)
ranef(model1)

### EXTRACT THE FITTED VALUES FROM THE BEST LMER MODEL ###
### YOU HAVE TO CHECK THE AIC CRITERIA ###
### FOR ALPHAT TRAIT, BEST MODEL IS 1    puede ser el modelo 1c  que un poquito menos 









ranef(model1)
rr1 = ranef(model1, condVar = TRUE)
plot1 = dotplot(rr1)[1]
plot2 = dotplot(rr1)[2]
plot3 = dotplot(rr1)[3]


print(VarCorr(model0),comp="Variance")
print(VarCorr(model1),comp="Variance")

fixef(model0)
fixef(model1)

ranef(model0)
ranef(model1)

### EXTRACT THE FITTED VALUES FROM THE BEST LMER MODEL ###
### YOU HAVE TO CHECK THE AIC CRITERIA ###
### FOR ALPHAT TRAIT, BEST MODEL IS 1    puede ser el modelo 1c  que un poquito menos 

fitted.FLOWERING = data.frame(fitted.values(model1))
colnames(fitted.FLOWERING)= c("fitted_FLOWERING")
head(fitted.FLOWERING)

hist(fitted.FLOWERING$fitted_FLOWERING)
hist(corr.data.FLOWERING$FLOWERING)

write.csv(fitted.FLOWERING,file ="fitted.FLOWERING_model1.csv" )



#### model for leaf####

model0L = lmer(LEAF ~ (1|RIL), data=corr.df)
model0bL = lmer(LEAF ~ RIL + (1|RIL), data = corr.df)

### add foto  (trataiento sd o ld)
#model1L = lmer(LEAF ~ (1|RIL) + (1|REP), data=corr.df)
#error de convergencia
#model1bL = lmer(LEAF ~ RIL + (1|foto), data=corr.df)
#model1cL = lmer(LEAF ~ RIL + foto + (1|RIL:foto), data=corr.df)

### ADD EFFECT of 'mesured_in' parameter###
model1L = lmer(LEAF ~ (1|RIL) + (1|REP), data=corr.df)
model1bL = lmer(LEAF ~ RIL + (1|REP), data=corr.df)
#model1c = lmer(LEAF ~ RIL + REP + (1|RIL:REP), data=corr.df)


# and its interaction with the genotype
#model2L = lmer(LEAF ~ (1|RIL) + (1|REP) + (1|RIL:REP), data=corr.df)
#model2bL = lmer(LEAF ~ RIL + (1|RIL) + (1|REP) + (1|RIL:REP), data=corr.df)
# errot de convergencia
### ADD THE EFFECT OF 'Treatment' parameter to model 1
#model3 = lmer(alphaT ~ (1|Genotype_name) + (1|Measured_in) + (1|Treatment), data=corr.df)
model3L = lmer(LEAF ~ (1|RIL) + (1|REP), data=corr.df)

# and its interaction with the genotype

#model4 = lmer(FLOWERING ~ (1|Genotype_name) + (1|Measured_in) + (1|Treatment) + (1|Genotype_name:Treatment), data=corr.df)

model4L = lmer(LEAF ~ (1|RIL) + (1|REP) + (1|foto), data=corr.df)

# add the interaction with the 'measured_in'
model4bL = lmer(LEAF ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto), data=corr.df)
model5L = lmer(LEAF ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto) + (1|RIL:REP), data=corr.df)
#model6 = lmer(FLOWERING ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto) + (1|RIL:REP) + (1|foto:REP), data=corr.df)
#no le gusta el modelo6 modelo falla en coverngncia 

#model7 = lmer(FLOWERING ~ (1|RIL)  + (1|RIL:foto) + (1|RIL:REP), data=corr.df)
#modelo sobre ajustado


model8L = lmer(LEAF ~ (1|RIL) + (1|foto) + (1|RIL:foto) + (1|RIL:REP), data=corr.df)


anova(model0L,model0bL,model1L,model1bL,model3L)

VarCorr(model0L)
VarCorr(model0bL)
VarCorr(model1L)
VarCorr(model1bL)
VarCorr(model3L)


VarCorr(model0L)
VarCorr(model0bL)
VarCorr(model1L)
VarCorr(model1bL)
VarCorr(model1cL)
VarCorr(model2L)
VarCorr(model3L)
VarCorr(model4L)
VarCorr(model4bL)
VarCorr(model5)
VarCorr(model8L)

ranef(model1bL)
rr1 = ranef(model1bL, condVar = TRUE)
plot1 = dotplot(rr1)[1]
plot2 = dotplot(rr1)[2]
plot3 = dotplot(rr1)[3]


print(VarCorr(model0),comp="Variance")
print(VarCorr(model1bL),comp="Variance")

fixef(model0)
fixef(model1bL)

ranef(model0)
ranef(model1L)

### EXTRACT THE FITTED VALUES FROM THE BEST LMER MODEL ###
### YOU HAVE TO CHECK THE AIC CRITERIA ###
### FOR ALPHAT TRAIT, BEST MODEL IS 1L

fitted.LEAF = data.frame(fitted.values(model1bL))
colnames(fitted.LEAF)= c("fitted_LEAF")
head(fitted.LEAF)

hist(fitted.LEAF$fitted_LEAF)
hist(corr.data.LEAF$LEAF)

write.csv(fitted.LEAF,file ="fitted.LEAF_model1.csv" )



#### model for petiole####
colnames(corr.df)
model0P = lmer(PETIOLE ~ (1|RIL), data=corr.df)
model0bP = lmer(PETIOLE ~ RIL + (1|RIL), data = corr.df)

### add foto  (trataiento sd o ld)
model1P = lmer(PETIOLE ~ (1|RIL) + (1|REP), data=corr.df)
model1bP = lmer(PETIOLE ~ RIL + (1|REP), data=corr.df)
#model1cP= lmer(PETIOLE ~ RIL + foto + (1|RIL:REP), data=corr.df)

### ADD EFFECT of 'mesured_in' parameter###
#model1 = lmer(FLOWERING ~ (1|RIL) + (1|REP), data=corr.df)
#model1b = lmer(FLOWERING ~ RIL + (1|REP), data=corr.df)
#model1c = lmer(FLOWERING ~ RIL + REP + (1|RIL:REP), data=corr.df)


# and its interaction with the genotype
model2P = lmer(PETIOLE ~ (1|RIL) + (1|REP) + (1|RIL:REP), data=corr.df)
model2bP = lmer(PETIOLE ~ RIL + (1|RIL) + (1|REP) + (1|RIL:REP), data=corr.df)

### ADD THE EFFECT OF 'Treatment' parameter to model 1
#model3 = lmer(alphaT ~ (1|Genotype_name) + (1|Measured_in) + (1|Treatment), data=corr.df)
model3P = lmer(PETIOLE ~ (1|RIL) + (1|REP), data=corr.df)

# and its interaction with the genotype

#model4 = lmer(FLOWERING ~ (1|Genotype_name) + (1|Measured_in) + (1|Treatment) + (1|Genotype_name:Treatment), data=corr.df)

model4P = lmer(PETIOLE ~ (1|RIL) + (1|REP) + (1|foto), data=corr.df)

# add the interaction with the 'measured_in'
model4bP = lmer(PETIOLE ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto), data=corr.df)
model5P = lmer(PETIOLE ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto) + (1|RIL:REP), data=corr.df)
#model6 = lmer(FLOWERING ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto) + (1|RIL:REP) + (1|foto:REP), data=corr.df)
#no le gusta el modelo6 modelo falla en coverngncia 

#model7 = lmer(FLOWERING ~ (1|RIL)  + (1|RIL:foto) + (1|RIL:REP), data=corr.df)
#modelo sobre ajustado


model8P = lmer(PETIOLE ~ (1|RIL) + (1|foto) + (1|RIL:foto) + (1|RIL:REP), data=corr.df)


anova(model0P,model0bP,model1P,model1bP,model3P)

VarCorr(model0P)
VarCorr(model0bP)
VarCorr(model1P)
VarCorr(model1bP)
VarCorr(model3P)



VarCorr(model0P)
VarCorr(model0bP)
VarCorr(model1P)
VarCorr(model1bP)
VarCorr(model1cP)
VarCorr(model2P)
VarCorr(model2bP)
VarCorr(model3P)
VarCorr(model4P)
VarCorr(model4bP)
VarCorr(model5P)
VarCorr(model8P)

ranef(model1bP)
rr1 = ranef(model1bP, condVar = TRUE)
plot1 = dotplot(rr1)[1]
plot2 = dotplot(rr1)[2]
plot3 = dotplot(rr1)[3]


print(VarCorr(modelb0),comp="Variance")
print(VarCorr(model1bP),comp="Variance")

fixef(model0bP)
fixef(model1bP)

ranef(model0bP)
ranef(model1bP)

### EXTRACT THE FITTED VALUES FROM THE BEST LMER MODEL ###
### YOU HAVE TO CHECK THE AIC CRITERIA ###
### FOR ALPHAT TRAIT, BEST MODEL IS 1P

fitted.PETIOLE= data.frame(fitted.values(model1bP))
colnames(fitted.PETIOLE)= c("fitted_petiole")
head(fitted.PETIOLE)
#View(fitted)

hist(fitted.PETIOLE$fitted_PETIOLE)
hist(corr.data.PETIOLE$PETIOLE)

write.csv(fitted.PETIOLE,file ="fitted.PETIOLE_model1.csv" )





#### MODEL for lenght####
colnames(corr.df)
model0LE = lmer(LENGHT~ (1|RIL), data=corr.df)
model0bLE = lmer(LENGHT ~ RIL + (1|RIL), data = corr.df)

### add foto  (trataiento sd o ld)
model1LE = lmer(LENGHT ~ (1|RIL) + (1|REP), data=corr.df)
model1bLE = lmer(LENGHT ~ RIL + (1|REP), data=corr.df)
#model1cLE = lmer(LENGHT~ RIL + foto + (1|RIL:REP), data=corr.df)

### ADD EFFECT of 'mesured_in' parameter###
#model1 = lmer(FLOWERING ~ (1|RIL) + (1|REP), data=corr.df)
#model1b = lmer(FLOWERING ~ RIL + (1|REP), data=corr.df)
#model1c = lmer(FLOWERING ~ RIL + REP + (1|RIL:REP), data=corr.df)


# and its interaction with the genotype
model2LE = lmer(LENGHT~ (1|RIL) + (1|foto) + (1|RIL:foto), data=corr.df)
model2bLE = lmer(LENGHT ~ RIL + (1|RIL) + (1|foto) + (1|RIL:foto), data=corr.df)

### ADD THE EFFECT OF 'Treatment' parameter to model 1
#model3 = lmer(alphaT ~ (1|Genotype_name) + (1|Measured_in) + (1|Treatment), data=corr.df)
model3LE = lmer( LENGHT ~ (1|RIL) + (1|REP), data=corr.df)

# and its interaction with the genotype

#model4 = lmer(FLOWERING ~ (1|Genotype_name) + (1|Measured_in) + (1|Treatment) + (1|Genotype_name:Treatment), data=corr.df)

model4LE = lmer( LENGHT ~ (1|RIL) + (1|REP) + (1|foto), data=corr.df)

# add the interaction with the 'measured_in'
model4bLE = lmer( LENGHT ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto), data=corr.df)
model5LE = lmer( LENGHT~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto) + (1|RIL:REP), data=corr.df)
#model6 = lmer(FLOWERING ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto) + (1|RIL:REP) + (1|foto:REP), data=corr.df)
#no le gusta el modelo6 modelo falla en coverngncia 

#model7 = lmer(FLOWERING ~ (1|RIL)  + (1|RIL:foto) + (1|RIL:REP), data=corr.df)
#modelo sobre ajustado


model8LE = lmer( LENGHT~ (1|RIL) + (1|foto) + (1|RIL:foto) + (1|RIL:REP), data=corr.df)


anova(model0LE,model0bLE,model1LE,model1bLE,model3LE)

VarCorr(model0LE)
VarCorr(model0bLE)
VarCorr(model1LE)
VarCorr(model1bLE)
VarCorr(model3LE)





VarCorr(model0LE)
VarCorr(model0bLE)
VarCorr(model1LE)
VarCorr(model1bLE)
VarCorr(model1cLE)
VarCorr(model2LE)
VarCorr(model2bLE)
VarCorr(model3LE)
VarCorr(model4LE)
VarCorr(model4bLE)
VarCorr(model5LE)
VarCorr(model8LE)

ranef(model1bLE)
rr1 = ranef(model1bLE, condVar = TRUE)
plot1 = dotplot(rr1)[1]
plot2 = dotplot(rr1)[2]
plot3 = dotplot(rr1)[3]


print(VarCorr(model0LE),comp="Variance")
print(VarCorr(model1bLE),comp="Variance")

fixef(model0LE)
fixef(model1bLE)

ranef(model0LE)
ranef(model1bLE)

### EXTRACT THE FITTED VALUES FROM THE BEST LMER MODEL ###
### YOU HAVE TO CHECK THE AIC CRITERIA ###
### FOR ALPHAT TRAIT, BEST MODEL IS 1


fitted.LENGHT= data.frame(fitted.values(model1bLE))
colnames(fitted.LENGHT)= c("fitted_LENGHT")
head(fitted.LENGHT)
#View(fitted)

hist(fitted.LENGHT$fitted_LENGHT)
hist(corr.data.LENGHT$LENGHT)

write.csv(fitted.LENGHT,file ="fitted.LENGHT_model1.csv" )



#### MODEL FOR WIGHT####

colnames(corr.df)
model0W = lmer(WIGHT ~ (1|RIL), data=corr.df)
model0bW = lmer(WIGHT ~ RIL + (1|RIL), data = corr.df)

### add foto  (trataiento sd o ld)
model1W = lmer(WIGHT ~ (1|RIL) + (1|REP), data=corr.df)
model1bW = lmer(WIGHT ~ RIL + (1|REP), data=corr.df)
model1cW = lmer(WIGHT ~ RIL + REP + (1|RIL:REP), data=corr.df)

### ADD EFFECT of 'mesured_in' parameter###
#model1 = lmer(FLOWERING ~ (1|RIL) + (1|REP), data=corr.df)
#model1b = lmer(FLOWERING ~ RIL + (1|REP), data=corr.df)
#model1c = lmer(FLOWERING ~ RIL + REP + (1|RIL:REP), data=corr.df)


# and its interaction with the genotype
model2W = lmer(WIGHT ~ (1|RIL) + (1|foto) + (1|RIL:foto), data=corr.df)
model2bW = lmer(WIGHT ~ RIL + (1|RIL) + (1|foto) + (1|RIL:foto), data=corr.df)

### ADD THE EFFECT OF 'Treatment' parameter to model 1
#model3 = lmer(alphaT ~ (1|Genotype_name) + (1|Measured_in) + (1|Treatment), data=corr.df)
model3W = lmer(WIGHT ~ (1|RIL) + (1|REP), data=corr.df)

# and its interaction with the genotype

#model4 = lmer(FLOWERING ~ (1|Genotype_name) + (1|Measured_in) + (1|Treatment) + (1|Genotype_name:Treatment), data=corr.df)

model4W = lmer(WIGHT ~ (1|RIL) + (1|REP) + (1|foto), data=corr.df)

# add the interaction with the 'measured_in'
model4bW = lmer(WIGHT ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto), data=corr.df)
model5W = lmer(WIGHT ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto) + (1|RIL:REP), data=corr.df)
#model6 = lmer(FLOWERING ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto) + (1|RIL:REP) + (1|foto:REP), data=corr.df)
#no le gusta el modelo6 modelo falla en coverngncia 

#model7 = lmer(FLOWERING ~ (1|RIL)  + (1|RIL:foto) + (1|RIL:REP), data=corr.df)
#modelo sobre ajustado


model8W = lmer(WIGHT ~ (1|RIL) + (1|foto) + (1|RIL:foto) + (1|RIL:REP), data=corr.df)


anova(model0W,model0bW,model1W,model1bW,model3W)

VarCorr(model0W)
VarCorr(model0bW)
VarCorr(model1W)
VarCorr(model1bW)
VarCorr(model3W)


VarCorr(model0W)
VarCorr(model0bW)
VarCorr(model1W)
VarCorr(model1bW)
VarCorr(model1cW)
VarCorr(model2W)
VarCorr(model2bW)
VarCorr(model3W)
VarCorr(model4W)
VarCorr(model4bW)
VarCorr(model5W)
VarCorr(model8W)

ranef(model1bW)
rr1 = ranef(model1bW, condVar = TRUE)
plot1 = dotplot(rr1)[1]
plot2 = dotplot(rr1)[2]
plot3 = dotplot(rr1)[3]


print(VarCorr(model0W),comp="Variance")
print(VarCorr(model1bW),comp="Variance")

fixef(model0W)
fixef(model1bW)

ranef(model0W)
ranef(model1bW)

fitted.WIGHT = data.frame(fitted.values(model1bW))
colnames(fitted.WIGHT)= c("fitted_WIGHT")
head(fitted.WIGHT)

hist(fitted.WIGHT$fitted_WIGHT)
hist(corr.data.WIGHT$WIGHT)

write.csv(fitted.WIGHT,file ="fitted.WIGHT_model1.csv" )
## all trait use model1+fisrt letter of trait names (col names)


#### MODEL FOR fr.osc####

colnames(corr.df)
model0W1 = lmer(fr.osc ~ (1|RIL), data=corr.df)
model0bW1 = lmer(fr.osc ~ RIL + (1|RIL), data = corr.df)

### add foto  (trataiento sd o ld)
model1W1 = lmer(fr.osc ~ (1|RIL) + (1|REP), data=corr.df)
model1bW1 = lmer(fr.osc ~ RIL + (1|REP), data=corr.df)
#model1cW = lmer(fr.osc ~ RIL + REP + (1|RIL:REP), data=corr.df)

### ADD EFFECT of 'mesured_in' parameter###
#model1 = lmer(FLOWERING ~ (1|RIL) + (1|REP), data=corr.df)
#model1b = lmer(FLOWERING ~ RIL + (1|REP), data=corr.df)
#model1c = lmer(FLOWERING ~ RIL + REP + (1|RIL:REP), data=corr.df)


# and its interaction with the genotype
model2W = lmer(WIGHT ~ (1|RIL) + (1|foto) + (1|RIL:foto), data=corr.df)
model2bW = lmer(WIGHT ~ RIL + (1|RIL) + (1|foto) + (1|RIL:foto), data=corr.df)

### ADD THE EFFECT OF 'Treatment' parameter to model 1
#model3 = lmer(alphaT ~ (1|Genotype_name) + (1|Measured_in) + (1|Treatment), data=corr.df)
model3W1 = lmer(fr.osc ~ (1|RIL) + (1|REP), data=corr.df)

# and its interaction with the genotype

#model4 = lmer(FLOWERING ~ (1|Genotype_name) + (1|Measured_in) + (1|Treatment) + (1|Genotype_name:Treatment), data=corr.df)

model4W = lmer(WIGHT ~ (1|RIL) + (1|REP) + (1|foto), data=corr.df)

# add the interaction with the 'measured_in'
model4bW = lmer(WIGHT ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto), data=corr.df)
model5W = lmer(WIGHT ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto) + (1|RIL:REP), data=corr.df)
#model6 = lmer(FLOWERING ~ (1|RIL) + (1|REP) + (1|foto) + (1|RIL:foto) + (1|RIL:REP) + (1|foto:REP), data=corr.df)
#no le gusta el modelo6 modelo falla en coverngncia 

#model7 = lmer(FLOWERING ~ (1|RIL)  + (1|RIL:foto) + (1|RIL:REP), data=corr.df)
#modelo sobre ajustado


model8W = lmer(WIGHT ~ (1|RIL) + (1|foto) + (1|RIL:foto) + (1|RIL:REP), data=corr.df)


anova(model0W1,model0bW1,model1W1,model1bW1,model3W1)

VarCorr(model0W1)
VarCorr(model0bW1)
VarCorr(model1W1)
VarCorr(model1bW1)
VarCorr(model3W1)


VarCorr(model0W)
VarCorr(model0bW)
VarCorr(model1W)
VarCorr(model1bW)
VarCorr(model1cW)
VarCorr(model2W)
VarCorr(model2bW)
VarCorr(model3W)
VarCorr(model4W)
VarCorr(model4bW)
VarCorr(model5W)
VarCorr(model8W)

ranef(model1bW1)
rr1 = ranef(model1bW1, condVar = TRUE)
plot1 = dotplot(rr1)[1]
plot2 = dotplot(rr1)[2]
plot3 = dotplot(rr1)[3]


print(VarCorr(model0W1),comp="Variance")
print(VarCorr(model1bW1),comp="Variance")

fixef(model0W1)
fixef(model1bW1)

ranef(model0W1)
ranef(model1bW1)

fitted.WIGHT = data.frame(fitted.values(model1bW1))
colnames(fitted.WIGHT)= c("fitted_WIGHT")
head(fitted.WIGHT)

hist(fitted.WIGHT$fitted_WIGHT)
hist(corr.data.WIGHT$WIGHT)

write.csv(fitted.WIGHT,file ="fitted.fr.osc_model1.csv" )




