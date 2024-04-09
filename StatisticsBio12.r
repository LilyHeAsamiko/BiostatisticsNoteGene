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
#!!!!!!!!!!!!!!!!!!Factor Analysis Exercises
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
library(Biobase)
library(GSE5859Subset)
data(GSE5859Subset)
#Suppose you want to make an MA plot of the first two samples y = geneExpression[,1:2]. Which of the following projections gives us the projection of y
#so that column2 versus column 1 is an MA plot?
  
#A.y(1/2–√12–√1/2–√−1/2–√)B.y(111−1)C.(111−1)yD.(111−1)y⊤
#A

#Say Y is M×N, in the SVD Y=UDV⊤ which of the following is not correct?
#A) DV⊤ are the new coordinates for the projection U⊤Y
#B) UD are the new coordinates for the projection YV
#C) D are the coordinates of the projection U⊤Y
#D) U⊤Y is a projection from an N-dimensional to M-dimensional subspace.
#C

#Define:
#y = geneExpression - rowMeans(geneExpression)
#Compute and plot an image of the correlation for each sample. Make two image plots of these correlations. In the first one, plot the correlation as image. In the second, order the samples by date and then plot an image of the correlation. The only difference in these plots is the order in which the samples are plotted.
library(Biobase)
#library(GSE5859)
#data(GSE5859)
load('GSE5859.rdata')
n <- nrow(pData(e))
o <- order(pData(e)$date)
Y=exprs(e)[,o]
cors=cor(Y-rowMeans(Y))
cols=colorRampPalette(brewer.pal(11,"RdBu"))(100)
par(mfrow = c(1,2))
image(1:n,1:n,cors,xaxt="n",yaxt="n",col=cols,xlab="",ylab="",zlim=c(-1,1),main='ordered by time')
image(1:n,1:n,cor(exprs(e)-rowMeans(exprs(e))),xaxt="n",yaxt="n",col=cols,xlab="",ylab="",zlim=c(-1,1),main='original correlations')
#Based on these plots, which of the following you would say is true?
#A) The samples appear to be completely independent of each other.
#B) Sex seems to be creating structures as evidenced by the two cluster of highly correlated samples.
#C) There appear to be only two factors completely driven by month.
#D) The fact that in the plot ordered by month we see two groups mainly driven by month, and within these we see subgroups driven by date, seems to suggest date more than month per se are the hidden factors.
#D
#Based on the correlation plots above, we could argue that there are at least two hidden factors. Using PCA estimate these two factors. Specifically, apply the svd to y and use the first two PCs as estimates.
#Which command gives us these estimates?
#A) pcs = svd(y)$v[1:2,]
#B) pcs = svd(y)$v[,1:2]
#C) pcs = svd(y)$u[,1:2]
#D) pcs = svd(y)$d[1:2]
#C
#Plot each of the estimated factors ordered by date. Use color to denote month. The first factor is clearly related to date. Which of the following appear to be most different according to this factor?
#A) June 23 and June 27
#B) Oct 07 and Oct 28
#C) June 10 and June 23
#D) June 15 and June 24
DMY = pData(e)$date
md = format(DMY,'%m%d')
D1 =which(md[o]=='0623')
D2 = which(md[o]=='0627')
D3 =which(md[o]=='1007')
D4 = which(md[o]=='1028')
D5 =which(md[o]=='0610')
D6 = which(md[o]=='0623')
D7 =which(md[o]=='0615')
D8 = which(md[o]=='0624')
y = Y-rowMeans(Y)
pcs = svd(y)$u

cols=colorRampPalette(brewer.pal(11,"RdBu"))(100)
par(mfrow = c(1,1))
image(1:n,1:n,cor(pc1-rowMeans(pcs)),xaxt="n",yaxt="n",col=cols,xlab="",ylab="",zlim=c(-1,1),main='pc1 correlations')
par(mfrow = c(2,2))
image(cor(pcs[,c(D1,D2)]-rowMeans(pcs[,c(D1,D2)])),main=paste0('mean of pc correlations',mean(cor(pcs[,c(D1,D2)]-rowMeans(pcs[,c(D1,D2)])))))
image(cor(pcs[,c(D3,D4)]-rowMeans(pcs[,c(D3,D4)])),main=paste0('mean of pc correlations',mean(cor(pcs[,c(D3,D4)]-rowMeans(pcs[,c(D3,D4)])))))
image(cor(pcs[,c(D5,D6)]-rowMeans(pcs[,c(D5,D6)])),main=paste0('mean of pc correlations',mean(cor(pcs[,c(D5,D6)]-rowMeans(pcs[,c(D5,D6)])))))
image(cor(pcs[,c(D7,D8)]-rowMeans(pcs[,c(D7,D8)])),main=paste0('mean of pc correlations',mean(cor(pcs[,c(D7,D8)]-rowMeans(pcs[,c(D7,D8)])))))
#A
#Use the svd function to obtain the principal components (PCs) for our detrended gene expression data y.
explainpco = svd(y)$d^2/sum(svd(y)$d^2)*100
explainpc = cumsum(svd(y)$d^2)/sum(svd(y)$d^2)*100
#How many PCs explain more than 10% of the variability?
par(mfrow = c(1,1))
plot(explainpco,main = 'percent explained by odered PC ')
length(which(explainpco>10))
#1
#Which PC most correlates (negative or positive correlation) with month?
mo = format(DMY,'%m')
om = order(DMY)
#What is this correlation (in absolute value)?
ym = Y[,om]-rowMeans(Y[,om])
mpcs = svd(ym)$u
mpcorder=order(abs(cor(mpcs-rowMeans(mpcs))))
which(mpcorder == length(mpcorder))
#43264
max(abs(cor(mpcs-rowMeans(mpcs))))
#1
#Which PC most correlates (negative or positive correlation) with sex?
#What is this correlation (in absolute value)?
#so = pData(e)$group
Y=geneExpression
so = sampleInfo$group
os = order(so)
ys = Y[,os]-rowMeans(Y[,os])
spcs = svd(ys)$u
spcorder=order(abs(cor(spcs-rowMeans(spcs))))
which(spcorder == length(spcorder))
#576
max(abs(cor(spcs-rowMeans(spcs))))
#1
#Now instead of using month, which we have shown does not quite describe the batch, add the two estimated factors s$v[,1:2] to the linear model we used above. Apply this model to each gene and compute q-values for the sex difference. How many q-values <
#0.1 for the sex comparison?
library(sva) #available from Bioconductor
mod <- model.matrix(~sex)
#cleandat <- ComBat(geneExpression,svd(ys)$v[,1:2],mod)
#Error in ComBat(geneExpression, svd(ys)$v[, 1:2], mod) : 
#  This version of ComBat only allows one batch variable
Y=geneExpression
mo = format(sampleInfo$date,'%m')
om = order(mo)
#What is this correlation (in absolute value)?
ym = Y[,om]-rowMeans(Y[,om])
cleandat <- ComBat(ym,svd(ym)$v[,1],mod)
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

#What proportion of the genes are on chromosomes X and Y?
X1 = model.matrix(~sex+svd(ym)$v[,1:2])
#How many of the q-values for the group comparison are now <0.1?
dtmp = Y[,om]
res1 <- t( sapply(1:nrow(dtmp),function(j){
  y <- dtmp[j,]
  fit1 <- lm(y~X1-1)
  summary(fit1)$coef[2,]
} ) )


##turn into data.frame so we can use the same code for plots as above
res1 <- data.frame(res1[,c(1,4)])
names(res1) <- c("dm","p.value")

par(mfrow=c(1,2))
hist(res1$p.value[which(!chr%in%c("chrX","chrY") )],main="",ylim=c(0,1300))

plot(res1$dm,-log10(res1$p.value),main = '~sex+svd$v[,1:2]')
points(res1$dm[which(chr=="chrX")],-log10(res0$p.value[which(chr=="chrX")]),col=1,pch=16)
points(res1$dm[which(chr=="chrY")],-log10(res0$p.value[which(chr=="chrY")]),col=2,pch=16,xlab="Effect size",ylab="-log10(p-value)")
legend("bottomright",c("chrX","chrY"),col=1:2,pch=16)
qvals1 <- qvalue(res1$p.value)$qvalue
index11 <- which(qvals1<0.1)
abline(h=-log10(max(res1$p.value[index11])))  
#Note correction.
#length(index11)/length(index)
#0.01470588
#With this new list, what propotion of chrX and chrY?
xpc1 = length(which(!is.na(match(index11,which(chr=="chrX")))))/length(index11)  
#0
ypc1 = length(which(!is.na(match(index11,which(chr=="chrY")))))/length(index11)  
#0.00325325
#How many on Y or X?
length(which(!is.na(match(index11,which(chr=="chrX")))))+length(which(!is.na(match(index11,which(chr=="chrY")))))
#188

#Adjusting with factor analysis exercises
#Exercises
#In this section we will use the sva function in the sva package (available from Bioconductor) and apply it to the following data:
library(sva)
library(Biobase)
library(GSE5859Subset)
data(GSE5859Subset)
#In a previous section we estimated factors using PCA, but we noted that the first factor was correlated with our outcome of interest:
s <- svd(geneExpression-rowMeans(geneExpression))
cor(sampleInfo$group,s$v[,1])
#0.6236858
#The svafit function estimates factors, but downweighs the genes that appear to correlate with the outcome of interest. It also tries to estimate the number of factors and returns the estimated factors like this:
sex = sampleInfo$group
mod = model.matrix(~sex)
svafit = sva(geneExpression,mod)
head(svafit$sv)
#Number of significant surrogate variables is:  5 
#Iteration (out of 5 ):1  2  3  4  5  > head(svafit$sv)
#[,1]        [,2]        [,3]       [,4]         [,5]
#[1,] -0.26862626 -0.03838109 -0.15306742  0.3007374 -0.210098159
#[2,] -0.06132157 -0.15755769  0.03538763 -0.0851655  0.063869257
#[3,] -0.12161818 -0.21766433  0.12624414 -0.2443445  0.099004174
#[4,]  0.30660574 -0.09657648  0.32034135  0.1680430 -0.620260643
#[5,] -0.01850853  0.18648507 -0.17931970 -0.4244993  0.007840835
#[6,]  0.36062840 -0.08542758  0.10726746 -0.1074114 -0.033204590
#The resulting estimated factors are not that different from the PCs.

for(i in 1:ncol(svafit$sv)){
  print( cor(s$v[,i],svafit$sv[,i]) )
}
#Now fit a linear model to each gene that instead of month includes these factors in the model. Use the qvalue function.
Xa = model.matrix(~sex+svafit$sv)
#How many of the q-values for the group comparison are now <0.1?
dtmp = Y[,om]
resa <- t( sapply(1:nrow(dtmp),function(j){
  y <- dtmp[j,]
  fit1 <- lm(y~Xa-1)
  summary(fit1)$coef[2,]
} ) )


##turn into data.frame so we can use the same code for plots as above
resa <- data.frame(resa[,c(1,4)])
names(resa) <- c("dm","p.value")

par(mfrow=c(1,2))
hist(resa$p.value[which(!chr%in%c("chrX","chrY") )],main="",ylim=c(0,1300))

plot(resa$dm,-log10(resa$p.value),main = '~sex+sva$sv')
points(resa$dm[which(chr=="chrX")],-log10(resa$p.value[which(chr=="chrX")]),col=1,pch=16)
points(resa$dm[which(chr=="chrY")],-log10(resa$p.value[which(chr=="chrY")]),col=2,pch=16,xlab="Effect size",ylab="-log10(p-value)")
legend("bottomright",c("chrX","chrY"),col=1:2,pch=16)
qvalsa <- qvalue(resa$p.value)$qvalue
indexa <- which(qvalsa<0.1)
abline(h=-log10(max(res1$p.value[indexa])))  
#Note correction.
#length(index11)/length(index)
#58.76471
#With this new list, what propotion of chrX and chrY?
xpca = length(which(!is.na(match(indexa,which(chr=="chrX")))))/length(indexa)  
#0.03941083
ypca = length(which(!is.na(match(indexa,which(chr=="chrY")))))/length(indexa)  
#0.001194268
#How many on Y or X?
length(which(!is.na(match(indexa,which(chr=="chrX")))))+length(which(!is.na(match(indexa,which(chr=="chrY")))))
#102

