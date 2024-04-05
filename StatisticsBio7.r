setwd("E:/")
library(dplyr)
#Hierarchical Models Exercises
#Exercises
#Load the following data (you can install it from Bioconductor) and extract the data matrix using exprs:
library(Biobase)
library(genefilter)
#library(SpikeInSubset)
#data(rma95)
#y <- exprs(rma95)
load("spikeInEDA.rda")
y<-int
#This dataset comes from an experiment in which RNA was obtained from the same background pool to create six replicate samples. Then RNA from 16 genes were artificially added in different quantities to each sample. These quantities (in picoMolars) and gene IDs are stored here:
#pData(rma95)
#These quantities where the same in the first three arrays and in the last three arrays. So we define two groups like this:
g <- factor(rep(0:1,each=3))
#and create an index of which rows are associated with the artificially added genes:
#spike <- rownames(y) %in% colnames(pData(rma95))
#Only these 16 genes are diferentially expressed since the six samples differ only due to sampling (they all come from the same background pool of RNA).
spike <- rownames(y) %in% colnames(spikeInDesign)
#Perform a t-test on each gene using the rowttest function.
tt00s={}
ttestr<-function(x){
  tt0<-t.test(x)
  array(list(tt0$p.value,tt0$statistic,tt0$conf.int))
}
for(r in 1:dim(spikeInDesign)[1]){tt00s[r]=list(c(t.test(spikeInDesign[r,])$p.value,t.test(spikeInDesign[r,])$statistic,t.test(spikeInDesign[r,])$conf.int))}
tt<-rowttests(int[,54:59],g)
tt0<-rowttests(t(spikeInDesign)[,54:59],g)
smallp <- with(tt, p.value < .01)
#What proportion of genes with a p-value < 0.01 (no multiple comparison correction) are not part of the artificially added (false positive)?
length(tt$p.value[which(tt$p.value<0.01)])/length(tt$p.value)
#0.625 [1]0.002205077
table(top50=rank(tt$p.value)<= 10, spike)
#spike
#top50    FALSE
#FALSE 201797
#TRUE      10
#Now compute the within group sample standard deviation for each gene (you can use group 1). Based on the p-value cut-off, split the genes into true positives, false positives, true negatives and false negatives. Create a boxplot comparing the sample SDs for each group. Which of the following best describes the boxplot?
#A) The standard deviation is similar across all groups.
#B) On average, the true negatives have much larger variability.
#C) The false negatives have larger variability.
#D) The false positives have smaller standard deviation.
#In the previous two questions, we observed results consistent with the fact that the random variability associated with the sample standard deviation leads to t-statistics that are large by chance.
boxplot(tt$statistic)
#The sample standard deviation we use in the t-test is an estimate and with just a pair of triplicate samples, the variability associated with the denominator in the t-test can be large.
#The following steps perform the basic limma analysis. We specify coef=2 because we are interested in the difference between groups, not the intercept. The eBayes step uses a hierarchical model that provides a new estimate of the gene specific standard error.
library(limma)
fit <- lmFit(y[,54:59], design=model.matrix(~ g))
colnames(coef(fit))
fit <- eBayes(fit)
#Here is a plot of the original, new, hierarchical models based estimate versus the sample based estimate:
sampleSD = fit$sigma
posteriorSD = sqrt(fit$s2.post)

plot(sampleSD,posteriorSD)  
abline(0,1)

#Which best describes what the hierarchical model does?
  
#A) Moves all the estimates of standard deviation closer to 0.12.
#B) Increases the estimates of standard deviation to increase t.
#C) Decreases the estimate of standard deviation.
#D) Decreases the effect size estimates.

#Use these new estimates of standard deviation in the denominator of the t-test and compute p-values. You can do it like this:
library(limma)
fit = lmFit(y[,54:59], design=model.matrix(~ g))
fit = eBayes(fit)
##second coefficient relates to diffences between group
pvals = fit$p.value[,2] 
#What proportion of genes with a p-value < 0.01 (no multiple comparison correction) are not part of the artificially added (false positive)?
#Compare to the previous volcano plot and notice that we no longer have small p-values for genes with small effect sizes.
length(pvals[which(pvals<0.01)])/length(pvals)
#0.002110928
cols <- ifelse(spike,"dodgerblue",ifelse(smallp,"red","black"))
par(mfrow=c(1,2))
with(tt, plot(-dm, -log10(p.value),cex=.8, pch=16,
              xlab="difference in means",
              col=cols))
abline(h=2,v=c(-.2,.2), lty=2)

limmares <- data.frame(dm=coef(fit)[,2], p.value=fit$p.value[,1])
with(limmares, plot(dm, -log10(p.value),cex=.8, pch=16,
                    col=cols,xlab="difference in means"))
abline(h=2,v=c(-.2,.2), lty=2)

##Distance exercises
#Exercises
#If you have not done so already, install the data package tissueGeneExpression:
library(devtools)
install_github("genomicsclass/tissuesGeneExpression")
#The data represents RNA expression levels for eight tissues, each with several biological replictes. We call samples that we consider to be from the same population, such as liver tissue from different individuals, biological replicates:
library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e)
head(tissue)
#How many biological replicates for hippocampus?
length(tissue[which(tissue=='hippocampus')]) 
#What is the distance between samples 3 and 45?
sqrt( crossprod(e[,3]-e[,45]) )
#What is the distance between gene 210486_at and 200805_at
sqrt( crossprod(e['210486_at',]-e['200805_at',]) )
#If I run the command (don’t run it!):
d = as.matrix( dist(e) )
#how many cells (number of rows times number of columns) will this matrix have?
  
#Compute the distance between all pair of samples:
d = dist( t(e) )
#Read the help file for dist.How many distances are stored in d? Hint: What is the length of d?
#[1] 17766
#Why is the answer to exercise 5 not ncol(e)^2?
#C) Because we take advantage of symmetry: only lower triangular matrix is stored thus only ncol(e)*(ncol(e)-1)/2 values.
ncol(e)*(ncol(e)-1)/2

library(rafalib)
library(MASS)
n <- 100
y <- t(mvrnorm(n,c(0,0), matrix(c(1,0.95,0.95,1),2,2)))
s <- svd(y)

#SVD exercises
#For these exercises we are again going to use:
library(tissuesGeneExpression)
data(tissuesGeneExpression)
#Before we start these exercises, it is important to reemphasize that when using the SVD, in practice the solution to SVD is not unique. This because UDV⊤=(−U)D(−V)⊤
#. In fact, we can flip the sign of each column of U and, as long as we also flip the respective column in V
#, we will arrive at the same solution. Here is an example:
s = svd(e)
signflips = sample(c(-1,1),ncol(e),replace=TRUE)
signflips
#Now we switch the sign of each column and check that we get the same answer. We do this using the function sweep. If x is a matrix and a is a vector then sweep(x,1,y,FUN="*") applies the function FUN to each row i FUN(x[i,],a[i]), in this case x[i,]*a[i]. If instead of 1 we use 2, sweep applies this to columns. To learn about sweep read ?sweep`.
newu= sweep(s$u,2,signflips,FUN="*")
newv= sweep(s$v,2,signflips,FUN="*" )
identical( s$u %*% diag(s$d) %*% t(s$v), newu %*% diag(s$d) %*% t(newv))
#This is important to know because different implementations of the SVD algorithm may give different signs, which can lead to the same code resulting in different answers when run in different computer systems.

#1.Compute the SVD of e
s = svd(e)
#Now compute the mean of each row:
m = rowMeans(e)
#What is the correlation between the first column of U and m?
cor(s$u[,1],m)  
#2In exercise 1, we saw how the first column relates to the mean of the rows of e. If we change these means, the distances between columns do not change. For example, changing the means does not change the distances:
newmeans = rnorm(nrow(e)) ##random values we will add to create new means
newe = e+newmeans ##we change the means
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(newe[,3]-newe[,45])) 
#So we might as well make the mean of each row 0, since it does not help us approximate the column distances. We will define y as the detrended e and recompute the SVD:
y = e - rowMeans(e)
s = svd(y)
#We showed that UDV⊤ is equal to y up to numerical error:
resid = y - s$u %*% diag(s$d) %*% t(s$v)
max(abs(resid))
#[1] 1.212808e-12
#The above can be made more efficient in two ways. First, using the crossprod and, second, not creating a diagonal matrix. In R, we can multiply matrices x by vector a. The result is a matrix with rows i equal to x[i,]*a[i]. Run the following example to see this.
x=matrix(rep(c(1,2),each=5),5,2)
x*c(1:5)
#which is equivalent to:
sweep(x,1,1:5,"*")
#This means that we don’t have to convert s$d into a matrix.
#Which of the following gives us the same as diag(s$d)%*%t(s$v) ?
#A) s$d %*% t(s$v)
#B) s$d * t(s$v)
#C) t(s$d * s$v)
#D) s$v * s$d
diag(s$d)%*%t(s$v)-s$d * t(s$v) 
#B
#3.If we define vd = t(s$d * t(s$v)), then which of the following is not the same UDV⊤:
#A) tcrossprod(s$u,vd)
#B) s$u %*% s$d * t(s$v)
#C) s$u %*% (s$d * t(s$v) )
#D) tcrossprod( t( s$d*t(s$u)) , s$v)
resid =tcrossprod( t( s$d*t(s$u)) , s$v) - s$u %*% diag(s$d) %*% t(s$v)
max(abs(resid))
#4.Let z = s$d * t(s$v). We showed derivation demonstrating that because U
#is orthogonal, the distance between e[,3] and e[,45] is the same as the distance between y[,3] and y[,45] . which is the same as vd[,3] and vd[,45]
z = s$d * t(s$v)
##d was deinfed in question 2.1.5
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(y[,3]-y[,45]))
sqrt(crossprod(z[,3]-z[,45]))
#Note that the columns z have 189 entries, compared to 22,215 for e.

#What is the difference, in absolute value, between the actual distance:
#and the approximation using only two dimensions of z ?
abs(sqrt(crossprod(e[,3]-e[,45])))-abs(sqrt(crossprod(y[,3]-y[,45])))
abs(sqrt(crossprod(e[,3]-e[,45])))-abs(sqrt(crossprod(z[,3]-z[,45])))


#How many dimensions do we need to use for the approximation in exercise 4 to be within 10%?
  
#Compute distances between sample 3 and all other samples.
dt3=rep(0,189)
for(col in c(1:2,4:189)){dt3[col]<-crossprod(e[,3]-e[,col])}

#Recompute this distance using the 2 dimensional approximation. What is the Spearman correlation between this approximate distance and the actual distance?
#The last exercise shows how just two dimensions can be useful to get a rough idea about the actual distances.
s=svd(e)
#resid = y - s$u %*% diag(s$d) %*% t(s$v)
y=s$u %*% diag(s$d) %*% t(s$v)
s3=rep(0,189)
for(col in 1:189){s3[col]<-crossprod(y[,3]-y[,col])}

z = s$d * t(s$v)
approxdt3=rep(0,189)
for(col in 1:189){approxdt3[col]<-crossprod(z[,3]-z[,col])}
#for(col in 1:189){approxdt3[col]<-crossprod(z[3,]-z[col,])}

spc=cor(approxdt3,s3,method='spearman')
spco=cor(approxdt3,dt3,method='spearman')
spcc=cor(s3,dt3,method='spearman')


#MDS exercises
#Using the z we computed in exercise 4 of the previous exercises:
library(tissuesGeneExpression)
data(tissuesGeneExpression)
y = e - rowMeans(e)
s = svd(y)
z = s$d * t(s$v)
#we can make an mds plot:
library(rafalib)
ftissue = factor(tissue)
par(mfrow=c(1,1))
plot(z[1,],z[2,],col=as.numeric(ftissue),main='tissuesGeneExpression')
legend("topleft",levels(ftissue),col=seq_along(ftissue),pch=1)
#Now run the function cmdscale on the original data:
d = dist(t(e))
mds = cmdscale(d)
#What is the absolute value of the correlation between the first dimension of z and the first dimension in mds?
abs(cor(z[1,],mds[,1]))
#[1] 1
#What is the absolute value of the correlation between the second dimension of z and the second dimension in mds?
abs(cor(z[2,],mds[,2]))
#[1] 1

#Load the following dataset:
library(GSE5859Subset)
data(GSE5859Subset)
#Compute the svd and compute z.
s = svd(geneExpression-rowMeans(geneExpression))
z = s$d * t(s$v)
#Which dimension of z most correlates with the outcome sampleInfo$group ?
corm= z[1,]
for (r in 1:dim(z)[1]){corm[r]=abs(cor(z[r,],sampleInfo$group))}
#What is this max correlation?
max(corm)  
#Which dimension of z has the second highest correlation with the outcome sampleInfo$group?
which(rank(corm)==dim(z)[1]-1)
#[1] 6
corm[2]
#[1] 0.1104236
#Note these measurements were made during two months:
sampleInfo$date
#We can extract the month this way:
month = format( sampleInfo$date, "%m")
mon = factor( month)
#Which dimension of z has the second highest correlation with the outcome month
cord= z[1,]
for (r in 1:dim(z)[1]){cord[r]=abs(cor(z[r,],as.numeric(mon)))}
which(rank(cord)==dim(z)[1]-1)
#[1] 2
#What is this correlation?
cord[2]
#[1] 0.4479168
#(Advanced) The same dimension is correlated with both the group and the date. The following are also correlated:
  
table(sampleInfo$g,month)

#So is this first dimension related directly to group or is it related only through the month? Note that the correlation with month is higher. This is related to batch effects which we will learn about later.
table(sampleInfo$group,mon)
#   mon
#   06 10
#0  9  3
#1  3  9
cord-corm
# [1]  0.2061057226  0.3374932128 -0.0209449786 -0.1681444202
#[5] -0.1981448995 -0.4529838412 -0.1082263448 -0.0969998949
#[9]  0.0737383996 -0.1355479943 -0.2031613911 -0.1184744707
#[13]  0.0429559807 -0.0927542517 -0.0523073205  0.0241160820
#[17] -0.0228196611  0.0002566111 -0.0108485914  0.0249934323
#[21] -0.0044775046  0.0513835084  0.0343106822 -0.1819532293
mean(cord)-mean(corm)
#[1] -0.0446848
#In exercise 3 we saw that one of the dimensions was highly correlated to the sampleInfo$group. Now take the 5th column of U
#and stratify by the gene chromosome. Remove chrUn and make a boxplot of the values of U5
#stratified by chromosome.
U5=s$u[,5]
fs=geneAnnotation$CHR
which(fs!='chrUn')
par(mfrow=c(2,1))
groups1 <- split(U5[which(fs!='chrUn')],geneExpression[which(fs!='chrUn'),5]) 
boxplot(groups1)
#Which chromosome looks different from the rest? Copy and paste the name as it appears in geneAnnotation.
ind={}
for (l in 1:length(groups1)){if(type(groups1[l])!='double[1]'){ind[l]<-groups1[l]}}
#Given the answer to the last exercise, any guesses as to what sampleInfo$group represents?
groups2 <- split(U5[which(fs!='chrUn')],round(geneExpression[which(fs!='chrUn'),5])) 
boxplot(groups2)
sampleInfo$group
#sd of expression month
library(openxlsx)
# 创建一个空白的工作簿
wb <- createWorkbook()

# 在工作簿中创建一个工作表
addWorksheet(wb, sheetName = "Data")

# 在工作表中写入数据
writeData(wb, sheet = "Data", x = spikeInDesign)

# 保存工作簿为Excel文件
saveWorkbook(wb, file = "spikeInDesign.xlsx")
saveWorkbook(wb, file = "spikeInDesign.csv")