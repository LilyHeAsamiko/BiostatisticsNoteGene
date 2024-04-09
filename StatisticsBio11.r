setRepositories(ind =c(1:7))

library(stats)
library(dplyr)
library('pvclust')
library('dendextend')
install.packages("matrixStats")
library(matrixStats)
if (!require("gplots")) install.packages("gplots")
library(gplots)
library('genefilter')
if(!require(caret)) install.packages('caret')
library(caret)  
if(!require(class)) install.packages('class')
library(class)  
#Confounding Exercises
#Exercises
#Load the admissions data from the dagdata package (which is available from the genomicsclass repository):
library(Biobase)
library(devtools)
#install_github("genomicsclass/dagdata")
setwd("E:/")
#library(dagdata)
load('admissions.rda')
#Adjusting with linear models exercises
#Exercises
#For the dataset we have been working with, models do not help due to the almost perfect confounding. This is one reason we created the subset dataset:
library(GSE5859Subset)
data(GSE5859Subset)
#Here we purposely confounded month and group (sex), but not completely:
sex = sampleInfo$group
month = factor( format(sampleInfo$date,"%m"))
table( sampleInfo$group, month)
chr = geneAnnotation$CHR
#Using the functions rowttests and qvalue compare the two groups. Because this is a smaller dataset which decreases our power, we will use the more lenient FDR cut-off of 10%.
library(qvalue)
res <- rowttests(geneExpression,as.factor( sampleInfo$group ))
par(mfrow = c(1,2))
hist(res$p.value[which(!chr%in%c("chrX","chrY") )],main="",ylim=c(0,1300))

plot(res$dm,-log10(res$p.value))
points(res$dm[which(chr=="chrX")],-log10(res$p.value[which(chr=="chrX")]),col=1,pch=16)
points(res$dm[which(chr=="chrY")],-log10(res$p.value[which(chr=="chrY")]),col=2,pch=16,xlab="Effect size",ylab="-log10(p-value)")
legend("bottomright",c("chrX","chrY"),col=1:2,pch=16)
qvals <- qvalue(res$p.value,fdr.level = 0.1)$qvalue
index <- which(qvals<0.1)
abline(h=-log10(max(res$p.value[index])))
#How many gene have q-values less than 0.1?
length(index)  
#59
#Note that sampleInfo$group here presents males and females. Thus, we expect differences to be in on chrY and, for genes that escape inactivation, chrX. We do not expect many autosomal genes to be different between males and females. This gives us an opportunity to evaluate false and true positives with experimental data. For example, we evaluate results using the proportion genes of the list that are on chrX or chrY.
XYno = which(chr%in%c("chrX","chrY") )
#For the list calculated above, what proportion of this list is on chrX or chrY?
length(XYno)/length(chr)  
#0.04787899
#We can also check how many of the chromosomes X and Y genes we detected as different. How many are on Y?
length(which(!is.na(match(index,which(chr=="chrY")))))  
#8
#Now for the autosomal genes (not on chrX and chrY) for which q-value < 0.1, perform a t-test comparing samples processed in June to those processed in October.
#index2 = which(!is.na(match(index,which(!chr%in%c("chrX","chrY")))))
#What proportion of these have p-values <0.05 ?
res2 <- rowttests(geneExpression[XYno,],as.factor(month))
#t.test(Pair(geneExpression[index2,which(month==month[1])],geneExpression[index2,which(month==month[2])]))  
#The above result shows that the great majority of the autosomal genes show differences due to processing data. This provides further evidence that confounding is resulting in false positives. So we are going to try to model the month effect to better estimate the sex effect. We are going to use a linear model:
#pttest<-function(dat) {
#  for (i in 1:dim(dat)[1]){
#  t.test(Pair(dat[i,month==month[1]],dat[i,month==month[2]]))
#}}
#ptts = pttest(geneExpression)
qres2p = qvalue(res2$p.value)
index2 = which(qres2p$qvalues<0.1)
g21 = geneExpression[index2,month == month[1]]
g22 = geneExpression[index2,month == month[2]]
ptts = rep(0,dim(geneExpression[index2,])[1])
for (i in 1:dim(geneExpression[index2,])[1]){
  ptts[i]=t.test(Pair(geneExpression[index2[i],month==month[1]],geneExpression[index2[i],month==month[2]]))$p.value
}
length(ptts<0.05)/length(index2)
length(which(!is.na(match(res2$p.value<0.05,index2))))/length(index2)
#0.7761194
#Which of the following creates the appropriate design matrix?
  
#A) X = model.matrix(~sex+ethnicity)
#B) X = cbind(sex,as.numeric(month))
#C) It can’t be done with one line.
#D) X = model.matrix(~sex+month)
#D
#Now use the X defined above, to fit a regression model using lm for each gene. You can obtain p-values for estimated parameters using summary. Here is an example

X = model.matrix(~sex+month)
i = 234
y = geneExpression[i,]
fit = lm(y~X)
summary(fit)$coef
#How many of the q-values for the group comparison are now <0.1?
res <- t( sapply(1:nrow(geneExpression),function(j){
  y <- geneExpression[j,]
  fit <- lm(y~X-1)
  summary(fit)$coef[2,]
} ) )


##turn into data.frame so we can use the same code for plots as above
res <- data.frame(res[,c(1,4)])
names(res) <- c("dm","p.value")

par(mfrow=c(1,2))
hist(res$p.value[which(!chr%in%c("chrX","chrY") )],main="",ylim=c(0,1300))

plot(res$dm,-log10(res$p.value),main = '~sex+month')
points(res$dm[which(chr=="chrX")],-log10(res$p.value[which(chr=="chrX")]),col=1,pch=16)
points(res$dm[which(chr=="chrY")],-log10(res$p.value[which(chr=="chrY")]),col=2,pch=16,xlab="Effect size",ylab="-log10(p-value)")
legend("bottomright",c("chrX","chrY"),col=1:2,pch=16)
qvals <- qvalue(res$p.value)$qvalue
index0 <- which(qvals<0.1)
abline(h=-log10(max(res$p.value[index0])))  
#Note the big drop from what we obtained without the correction.
length(index0)/length(index)
#0.2881356
#With this new list, what proportion of these are chrX and chrY?
xpc = length(which(!is.na(match(index0,which(chr=="chrX")))))/length(index0)  
#0.52941
ypc = length(which(!is.na(match(index0,which(chr=="chrY")))))/length(index0)  
#0.35294#Notice the big improvement.
#How many on Y or X?
length(which(!is.na(match(index0,which(chr=="chrX")))))+length(which(!is.na(match(index0,which(chr=="chrY")))))
#14
#Now from the linear model above, extract the p-values related to the coefficient representing the October versus June differences using the same linear model.
#How many of the q-values for the month comparison are now <0.1?

X0 = model.matrix(~sex)
i = 234
y = geneExpression[i,]
fit = lm(y~X)
summary(fit)$coef
#How many of the q-values for the group comparison are now <0.1?
res0 <- t( sapply(1:nrow(geneExpression),function(j){
  y <- geneExpression[j,]
  fit0 <- lm(y~X0-1)
  summary(fit0)$coef[2,]
} ) )


##turn into data.frame so we can use the same code for plots as above
res0 <- data.frame(res0[,c(1,4)])
names(res0) <- c("dm","p.value")

par(mfrow=c(1,2))
hist(res0$p.value[which(!chr%in%c("chrX","chrY") )],main="",ylim=c(0,1300))

plot(res0$dm,-log10(res0$p.value),main = '~sex')
points(res0$dm[which(chr=="chrX")],-log10(res0$p.value[which(chr=="chrX")]),col=1,pch=16)
points(res0$dm[which(chr=="chrY")],-log10(res0$p.value[which(chr=="chrY")]),col=2,pch=16,xlab="Effect size",ylab="-log10(p-value)")
legend("bottomright",c("chrX","chrY"),col=1:2,pch=16)
qvals0 <- qvalue(res0$p.value)$qvalue
index00 <- which(qvals0<0.1)
abline(h=-log10(max(res0$p.value[index00])))  
#Note correction.
length(index00)/length(index)
#1
#With this new list, what propotion of chrX and chrY?
xpc0 = length(which(!is.na(match(index00,which(chr=="chrX")))))/length(index00)  
#0.20338983
ypc0 = length(which(!is.na(match(index00,which(chr=="chrY")))))/length(index00)  
#0.13559322
#How many on Y or X?
length(which(!is.na(match(index00,which(chr=="chrX")))))+length(which(!is.na(match(index00,which(chr=="chrY")))))
#20

#This approach is basically the approach implemented by Combat.
library(sva) #available from Bioconductor
mod <- model.matrix(~sex)
cleandat <- ComBat(geneExpression,month,mod)
## Found 2 batches
## Adjusting for 1 covariate(s) or covariate level(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
#Then the results can be used to fit a model with our variable of interest:
res<-genefilter::rowttests(cleandat,factor(sex))
#In this case, the results are less specific than what we obtain by fitting the simple linear model:
par(mfrow = c(1,2))
hist(res$p.value[which(!chr%in%c("chrX","chrY") )],main="~sex,usingcombat",ylim=c(0,1300))

plot(res$dm,-log10(res$p.value))
points(res$dm[which(chr=="chrX")],-log10(res$p.value[which(chr=="chrX")]),col=1,pch=16)
points(res$dm[which(chr=="chrY")],-log10(res$p.value[which(chr=="chrY")]),col=2,pch=16,xlab="Effect size",ylab="-log10(p-value)")
legend("bottomright",c("chrX","chrY"),col=1:2,pch=16)
qvals <- qvalue(res$p.value)$qvalue
index <- which(qvals<0.1)
abline(h=-log10(max(res$p.value[index])))

cat("Total genes with q-value < 0.1: ",length(index),"\n",
    "Number of selected genes on chrY: ", sum(chr[index]=="chrY",na.rm=TRUE),"\n",
    "Number of selected genes on chrX: ", sum(chr[index]=="chrX",na.rm=TRUE),sep="")
## Total genes with q-value < 0.1: 68
## Number of selected genes on chrY: 8
## Number of selected genes on chrX: 16
