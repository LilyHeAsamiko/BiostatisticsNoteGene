library(stats)
library(dplyr)
set.seed(1)
m = 10000
n = 24
x = matrix(rnorm(m*n),m,n)
colnames(x)=1:n
hc1 = hclust(dist(x),method='ward')
plot(hc1)

library('pvclust')
library('dendextend')
m = 10000
n = 24
x = matrix(rnorm(m*n),m,n)
colnames(x)=1:n
hc = pvclust(x,method.hclust='ward')
par(mfrow=c(1,2))
plot(hc)
ct = cutree(as.dendrogram(hc1),h=143)
hist(ct)
#cutree(as.dendrogram(hc),h=143)
#monte carlo
sd(rnorm(length(ct),mean(ct),sd(ct)))
sd(runif(length(ct),min(ct),max(ct)))

#Run kmeans with 4 centers for the blood RNA data:
  
library(GSE5859Subset)
data(GSE5859Subset)
plot(as.dendrogram(hc1 = hclust(dist(t(geneExpression)))))
month = format( sampleInfo$date, "%m")
mon = factor( month)
year = format( sampleInfo$date, "%y")
yea = factor( year)
#Set the seed to 10, set.seed(10) right before running kmeans with 5 centers.
set.seed(10)
#Explore the relationship of clusters and information in sampleInfo. Which of the following best describes what you find?
ct = cutree(hc1,k=5)
#plot(as.dendrogram(ct))
table(ct,sampleInfo$group)
#ct  0 1
#  1 8 3
#  2 0 5
#  3 2 3
#  4 1 1
#  5 1 0
table(sampleInfo$group,mon)
#mon
#  06 10
#0  9  3
#1  3  9
table(sampleInfo$group,sampleInfo$ethnicity)
#  ASN CEU HAN
#0  11   1   0
#1  12   0   0
#D

#A) sampleInfo$group is driving the clusters as the 0s and 1s are in completely different clusters.
#B) The year is driving the clusters.
#C) Date is driving the clusters.
#D) The clusters don’t depend on any of the column of sampleInfo

#Load the data:
  
library(GSE5859Subset)
data(GSE5859Subset)
#Pick the 25 genes with the highest across sample variance. This function might help:
  
install.packages("matrixStats")
library(matrixStats)
?rowMads ##we use mads due to a outlier sample
#Use heatmap.2 to make a heatmap showing the sampleInfo$group with color, the date as labels, the rows labelled with chromosome, and scaling the rows.
rmds = rowMads(matrixStats)
plot(rmds)
hist(rmds)
#indx = filter(rank(rmds),length(mds)-25+1:length(mds)) %>% unlist
#indx = rank(rmds)[which(rank(rmds) in length(rmds)-25+1:length(rmds))]
geneAnnotation['Rrmds'] = rank(rmds)#rank with increasing order
gene25 = filter(geneAnnotation, Rrmds %in% c(length(rmds)-25+1:length(rmds)))
addData = data.frame(geneExpression[gene25$PROBEID,])
addData['PROBEID'] = gene25$PROBEID
info25 = cbind(sampleInfo,t(data.frame(geneExpression[gene25$PROBEID,])))
geneExpression25 = merge(gene25,addData,by = 'PROBEID')

geneExpression25['DATE'] = sampleInfo$date
#xlab strCol = ,ylab strRow
if (!require("gplots")) install.packages("gplots")
library(gplots)
#data <- read.table(".txt",header=TRUE,row.names = 1)
#data <- as.matrix(geneExpression25[,6:29])
data <- info25[,4:29] #24*26
#data.rownames = info25$date
#data.colnames = c('group',geneExpression25$CHR)
rownames(data) = paste0(info25$filename,info25$date)
colnames(data) = c('group',geneExpression25$CHR)

coul <- colorRampPalette(colors = c("#4A73EE","white","#F46161"))(100)

heatmap.2(as.matrix(data),scale = "row",col = coul,
          main = 'hotplot',xlab = 'Chromosomes', ylab = 'Date',
          margins = c(10,20),  #行名和列名的边距。用于调整热图大小
          ColSideColors = rainbow(ncol(data)),  #注释列的水平边栏的颜色向量（行也类似）。
          
          Rowv = FALSE,   #determines if and how the row dendrogram should be reordered.即是否对行聚类并重新排序。取值为TRUE、FALSE或整数向量。
          Colv = FALSE,   #determines if and how the column dendrogram should be reordered. 是否对列进行聚类并重新排序。取值为TRUE、FALSE或“Rowv”（表示与行的处理方式相同）
          dendrogram="none", #whether to draw 'none', 'row', 'column' or 'both' dendrograms.即是否画出聚类树
          
          colCol = c(rep('red',ncol(data))), #列标签的颜色向量（行也类似）
          
#          rowsep = c(2,4),colsep = 4, #在指定的行/列插入gap
#          sepcolor = "white",  #gap的颜色
#          sepwidth=c(0.1,0.1),  #行、列中gap的大小
          
          trace = 'column',  # 是否在c("column","row","both","none")绘制实线,线越靠近单元格四侧，值越大或小。
          tracecol="green",  #实线颜色
          linecol='blue'   #中间虚线颜色
)

#What do we learn from this heatmap?
  
#A) The data appears as if it was generated by rnorm.
#B) Some genes in chr1 are very variable.
#C) A group of chrY genes are higher in group 0 and appear to drive the clustering. Within those clusters there appears to be clustering by month.
#D) A group of chrY genes are higher in October compared to June and appear to drive the clustering. Within those clusters there appears to be clustering by samplInfo$group.

#Create a large data set of random data that is completely independent of sampleInfo$group like this:
library('genefilter')

set.seed(17)
m = nrow(geneExpression)
n = ncol(geneExpression)
x = matrix(rnorm(m*n),m,n)
g = factor(sampleInfo$g )
#Create two heatmaps with these data. Show the group g either with labels or colors. First, take the 50 genes with smallest p-values obtained with rowttests. Then, take the 50 genes with largest standard deviations.
rtts = rowttests(geneExpression,g)
Rrtts = rank(rtts$p.value)
Rrttssd = rank(rtts$statistic)
info50 = cbind(rtts,Rrtts,Rrttssd)
#gene501 = geneExpression[filter(info50,Rrtts %in% c(1:50)),]
#gene502 = geneExpression[filter(info50,Rrttssd %in% c(length(Rrttssd)-49:length(Rrttssd))),]
gene501 = geneExpression[rownames(filter(info50,Rrtts %in% c(1:50))),]
gene502 = geneExpression[rownames(filter(info50,Rrttssd %in% c((length(info50$Rrttssd)-49):length(info50$Rrttssd)))),]
#par(mfrow = c(1,2))
coul <- colorRampPalette(colors = c("#4A73EE","white","#F46161"))(100)

heatmap.2(as.matrix(gene501),scale = "row",col = coul,
          main = 'hotplot',xlab = 'Person', ylab = 'Probeid',
          margins = c(10,20),  #行名和列名的边距。用于调整热图大小
          ColSideColors = rainbow(ncol(gene501)),  #注释列的水平边栏的颜色向量（行也类似）。
          
          Rowv = FALSE,   #determines if and how the row dendrogram should be reordered.即是否对行聚类并重新排序。取值为TRUE、FALSE或整数向量。
          Colv = FALSE,   #determines if and how the column dendrogram should be reordered. 是否对列进行聚类并重新排序。取值为TRUE、FALSE或“Rowv”（表示与行的处理方式相同）
          dendrogram="none", #whether to draw 'none', 'row', 'column' or 'both' dendrograms.即是否画出聚类树
          
          colCol = c(rep('red',ncol(gene501))), #列标签的颜色向量（行也类似）
          
          #          rowsep = c(2,4),colsep = 4, #在指定的行/列插入gap
          #          sepcolor = "white",  #gap的颜色
          #          sepwidth=c(0.1,0.1),  #行、列中gap的大小
          
          trace = 'column',  # 是否在c("column","row","both","none")绘制实线,线越靠近单元格四侧，值越大或小。
          tracecol="green",  #实线颜色
          linecol='blue'   #中间虚线颜色
)

coul <- colorRampPalette(colors = c("#4A73EE","white","#F46161"))(100)

heatmap.2(as.matrix(gene502),scale = "row",col = coul,
          main = 'hotplot',xlab = 'Person', ylab = 'Probeid',
          margins = c(10,20),  #行名和列名的边距。用于调整热图大小
          ColSideColors = rainbow(ncol(gene502)),  #注释列的水平边栏的颜色向量（行也类似）。
          
          Rowv = FALSE,   #determines if and how the row dendrogram should be reordered.即是否对行聚类并重新排序。取值为TRUE、FALSE或整数向量。
          Colv = FALSE,   #determines if and how the column dendrogram should be reordered. 是否对列进行聚类并重新排序。取值为TRUE、FALSE或“Rowv”（表示与行的处理方式相同）
          dendrogram="none", #whether to draw 'none', 'row', 'column' or 'both' dendrograms.即是否画出聚类树
          
          colCol = c(rep('red',ncol(gene502))), #列标签的颜色向量（行也类似）
          
          #          rowsep = c(2,4),colsep = 4, #在指定的行/列插入gap
          #          sepcolor = "white",  #gap的颜色
          #          sepwidth=c(0.1,0.1),  #行、列中gap的大小
          
          trace = 'column',  # 是否在c("column","row","both","none")绘制实线,线越靠近单元格四侧，值越大或小。
          tracecol="green",  #实线颜色
          linecol='blue'   #中间虚线颜色
)


#Which of the following statements is true?
  
#A) There is no relationship between g and x, but with 8,793 tests some will appear significant by chance. Selecting genes with the t-test gives us a deceiving result.
#B) These two techniques produced similar heatmaps.
#C) Selecting genes with the t-test is a better technique since it permits us to detect the two groups. It appears to find hidden signals.
#D) The genes with the largest standard deviation add variability to the plot and do not let us find the differences between the two groups.

#C

#Conditional Expectations Exercises
#Exercises
#Throughout these exercises it will be useful to remember that when our data are 0s and 1s, probabilities and expectations are the same thing. We can do the math, but here is some R code:
n = 1000
y = rbinom(n,1,0.25)
##proportion of ones Pr(Y)
sum(y==1)/length(y)
##expectaion of Y
mean(y)
#Generate some random data to imitate heights for men (0) and women (1):
n = 10000
set.seed(1)
men = rnorm(n,176,7) #height in centimeters
women = rnorm(n,162,7) #height in centimeters
y = c(rep(0,n),rep(1,n))#men woman labels
x = round(c(men,women))
##mix it up
ind = sample(seq(along=y))
y = y[ind]
x = x[ind]
#Using the data generated above, what is the $$E(Y	X=176)$?
fit = lm(y~x)
lb = fit$coefficient[1]+fit$coefficient[2]*176
#0.2494172
round(lb)
#0 
#Now make a plot of $$E(Y	X=x)$$ for x=seq(160,178) using the data generated in exercise 1.
#If you are predicting female or male based on height and want your probability of success to be larger than 0.5, what is the largest height where you predict female ?
par(mfrow=c(1,2))
plot(x,y,xlab='mixedHeight',ylab='mixedLabel',main=paste0('correlation=',signif(cor(x,y),2)))
abline(v=c(160,178),col ='red')
abline(fit,col = 'green')

x0=seq(160,178) 
hist(round(y[x==x0]),xlab='Label of mixedHeights',main = '',xlim =range(y))
abline(v = fit$coefficient[1]+fit$coefficient[2]*176)

filter(xx, )
xx = rank(x)
fid = rep(0,20000) #E(X=?|Y=1)
xn = 1
fidx = sum(y[which(x==min(x))])/length(y[which(x==min(x))])
#indx = which(xx)
#xtmp = xx[a]
#fid[a] = (round(y[xtmp])==1)-sum(y[which(x %in% xtmp)]==1))/sum(y[which(x %in% xtmp)]==1)
#fid[1] = 
for (xtmp in c((min(x)+1):max(x))){
  if (length(y[which(x==xtmp)])){
    fidx = c(fidx,sum(y[which(x==xtmp)])/length(y[which(x==xtmp)]))
    xn = c(xn,length(y[which(x==xtmp)]))}
  else{
    fidx = c(fidx,0)
    xn = c(xn,0)}
}
#
plot(xn,fidx) 
plot(fidx)
which(fidx>0.5)
#[1]  1  5  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28
#[25] 29 30 31 32 33 34 35 36 37
xn[which(fidx>0.5)]
# [1]   1   1   3   1   3   4   5  19  19  22  34  57  60 114 116 164 209 258
#[19] 319 377 393 498 532 535 610 647 644 675 674 700 664 700 686
min(x)+max(which(fidx>0.5))-1
#168

#Smoothing exercises
#Exercises
#Generate the following data:
n = 10000
set.seed(1)
men = rnorm(n,176,7) #height in centimeters
women = rnorm(n,162,7) #height in centimeters
y = c(rep(0,n),rep(1,n))
x = round(c(men,women))
##mix it up
ind = sample(seq(along=y))
y = y[ind]
x = x[ind]
#Set the seed at 5, set.seed(5) and take a random sample of 250 from:
set.seed(5)
N = 250
ind = sample(length(y),N)
Y = y[ind]
X = x[ind]
#Use loess to estimate $$f(x)=E(Y	X=x)usingthedefaultparameters.Whatisthepredicted
#f(168)$?
#The loess estimate above is a random variable. We can compute standard errors for it. Here we use Monte Carlo to demonstrate that it is a random variable. Use Monte Carlo simulation to estimate the standard error of your estimate of f(168).
Smoothing exercises
Exercises
Generate the following data:
  
  n = 10000
set.seed(1)
men = rnorm(n,176,7) #height in centimeters
women = rnorm(n,162,7) #height in centimeters
y = c(rep(0,n),rep(1,n))
x = round(c(men,women))
##mix it up
ind = sample(seq(along=y))
y = y[ind]
x = x[ind]
Set the seed at 5, set.seed(5) and take a random sample of 250 from:
  
  set.seed(5)
N = 250
ind = sample(length(y),N)
Y = y[ind]
X = x[ind]
Use loess to estimate $$f(x)=E(Y	X=x)usingthedefaultparameters.Whatisthepredicted
f(168)$?
  The loess estimate above is a random variable. We can compute standard errors for it. Here we use Monte Carlo to demonstrate that it is a random variable. Use Monte Carlo simulation to estimate the standard error of your estimate of f(168)
.

Set the seed to 5, set.seed(5) and perform 10000 simulations and report the SE of the loess based estimate.

Smoothing exercises
Exercises
Generate the following data:
  
  n = 10000
set.seed(1)
men = rnorm(n,176,7) #height in centimeters
women = rnorm(n,162,7) #height in centimeters
y = c(rep(0,n),rep(1,n))
x = round(c(men,women))
##mix it up
ind = sample(seq(along=y))
y = y[ind]
x = x[ind]
Set the seed at 5, set.seed(5) and take a random sample of 250 from:
  
  set.seed(5)
N = 250
ind = sample(length(y),N)
Y = y[ind]
X = x[ind]
Use loess to estimate $$f(x)=E(Y	X=x)usingthedefaultparameters.Whatisthepredicted
f(168)$?
  The loess estimate above is a random variable. We can compute standard errors for it. Here we use Monte Carlo to demonstrate that it is a random variable. Use Monte Carlo simulation to estimate the standard error of your estimate of f(168)
.

Set the seed to 5, set.seed(5) and perform 10000 simulations and report the SE of the loess based estimate.

Smoothing exercises
Exercises
Generate the following data:
  
  n = 10000
set.seed(1)
men = rnorm(n,176,7) #height in centimeters
women = rnorm(n,162,7) #height in centimeters
y = c(rep(0,n),rep(1,n))
x = round(c(men,women))
##mix it up
ind = sample(seq(along=y))
y = y[ind]
x = x[ind]
Set the seed at 5, set.seed(5) and take a random sample of 250 from:
  
  set.seed(5)
N = 250
ind = sample(length(y),N)
Y = y[ind]
X = x[ind]
Use loess to estimate $$f(x)=E(Y	X=x)usingthedefaultparameters.Whatisthepredicted
f(168)$?
  The loess estimate above is a random variable. We can compute standard errors for it. Here we use Monte Carlo to demonstrate that it is a random variable. Use Monte Carlo simulation to estimate the standard error of your estimate of f(168)
.

Set the seed to 5, set.seed(5) and perform 10000 simulations and report the SE of the loess based estimate.

gHKJ %>% %>% %>% 
#Set the seed to 5, set.seed(5) and perform 10000 simulations and report the SE of the loess based estimate.

