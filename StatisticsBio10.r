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
#data(admissions)
#Familiarize yourself with this table:
admissions
#Let’s compute the proportion of men who were accepted:
index = which(admissions$Gender==1)
accepted= sum(admissions$Number[index] * admissions$Percent[index]/100)
noaccepted = sum(admissions$Number[index] * (1-admissions$Percent[index]/100))
applied = sum(admissions$Number[index])
accepted/applied
#0.445191
#What is the proportion of women that were accepted?
index0 = which(admissions$Gender==0)
accepted0= sum(admissions$Number[index0] * admissions$Percent[index0]/100)
applied0 = sum(admissions$Number[index0])
noaccepted0 = sum(admissions$Number[index0] * (1-admissions$Percent[index0]/100))
accepted0/applied0
#0.3033351
#Now that we have observed different acceptance rates between genders, test for the significance of this result.
matrix(c(accepted,noaccepted,accepted0,noaccepted0),2,2)
#If you perform an independence test, what is the p-value?
chisq.test(matrix(c(accepted,noaccepted,accepted0,noaccepted0),2,2))  
#Pearson's Chi-squared test with Yates' continuity correction
#data:  matrix(c(accepted, noaccepted, accepted0, noaccepted0), 2, 2)
#X-squared = 91.895, df = 1, p-value < 2.2e-16
#This difference actually led to a lawsuit.

#Now notice that looking at the data by major, the differences disappear.
admissions
#How can this be? This is referred to as Simpson’s Paradox. In the following questions we will try to decipher why this is happening.

#We can quantify how “hard” a major is by using the percent of students that were accepted. Compute the percent that were accepted (regardless of gender) to each major and call this vector H.
A = which(admissions$Major=='A')
HA= sum(admissions$Number[A] * admissions$Percent[A]/100)
B = which(admissions$Major=='B')
HB= sum(admissions$Number[B] * admissions$Percent[B]/100)
C = which(admissions$Major=='C')
HC= sum(admissions$Number[C] * admissions$Percent[C]/100)
D = which(admissions$Major=='D')
HD= sum(admissions$Number[D] * admissions$Percent[D]/100)
E = which(admissions$Major=='E')
HE= sum(admissions$Number[E] * admissions$Percent[E]/100)
FF = which(admissions$Major=='F')
HF= sum(admissions$Number[FF] * admissions$Percent[FF]/100)
hardt =data.frame(HA,HB,HC,HD,HE,HF)
#Which is the hardest major?
colnames(hardt)[which(hardt == min(hardt))] 
#HF
#What proportion is accepted for this major?
min(hardt)
#46.25
#For men, what is the correlation between the number of applications across majors and H?
MA = which(admissions$Major=='A' & admissions$Gender==1)
MHA= sum(admissions$Number[MA] * admissions$Percent[MA]/100)
MB = which(admissions$Major=='B' & admissions$Gender==1)
MHB= sum(admissions$Number[MB] * admissions$Percent[MB]/100)
MC = which(admissions$Major=='C' & admissions$Gender==1)
MHC= sum(admissions$Number[MC] * admissions$Percent[MC]/100)
MD = which(admissions$Major=='D' & admissions$Gender==1)
MHD= sum(admissions$Number[MD] * admissions$Percent[MD]/100)
ME = which(admissions$Major=='E' & admissions$Gender==1)
MHE= sum(admissions$Number[ME] * admissions$Percent[ME]/100)
MF = which(admissions$Major=='F' & admissions$Gender==1)
MHF= sum(admissions$Number[MF] * admissions$Percent[MF]/100)

#For women, what is the correlation between the number of applications across majors and H?
FA = which(admissions$Major=='A' & admissions$Gender==0)
FHA= sum(admissions$Number[FA] * admissions$Percent[FA]/100)
FB = which(admissions$Major=='B' & admissions$Gender==0)
FHB= sum(admissions$Number[FB] * admissions$Percent[FB]/100)
FC = which(admissions$Major=='C' & admissions$Gender==0)
FHC= sum(admissions$Number[FC] * admissions$Percent[FC]/100)
FD = which(admissions$Major=='D' & admissions$Gender==0)
FHD= sum(admissions$Number[FD] * admissions$Percent[FD]/100)
FE = which(admissions$Major=='E' & admissions$Gender==0)
FHE= sum(admissions$Number[FE] * admissions$Percent[FE]/100)
FFF = which(admissions$Major=='F' & admissions$Gender==0)
FHF= sum(admissions$Number[FFF] * admissions$Percent[FF]/100)

#Given the answers to the above, which best explains the differences in admission percentages when we combine majors?
  
#A) We made a coding mistake when computing the overall admissions percentages.
#B) There were more total number of women applications which made the denominator much bigger.
#C) There is confounding between gender and preference for “hard” majors: females are more likely to apply to harder majors.
#D) The sample size for the individual majors was not large enough to draw the correct conclusion.

#EDA with PCA exercises
#Exercises
#We will use the Bioconductor package Biobase which you can install with install_bioc function from rafalib:
#Load the data for this gene expression dataset:
  
library(Biobase)
load('GSE5859.rda')
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)
#suppressMessages(library(GEOquery))
#library(Biobase)
#options( 'download.file.method.GEOquery' = 'libcurl' )
#options(timeout = 200000)
#gds <- getGEO("GSE86505" , destdir = "./",AnnotGPL = FALSE,getGPL = FALSE)#输出路径为默认，不下载大的GPL数据
#save(gset,file = 'GSE86505.gset.Rdata')
#install.packages('AnnoProbe')
#library(AnnoProbe)
#更新镜像库
#devtools::install_git("https://gitee.com/jmzeng/GEOmirror")
#使用中国镜像下载GEO数据
#gset <- AnnoProbe::geoChina(gse='GSE87211', mirror = 'tencent', destdir = '.')
#此处mirror仅有企鹅源
#devtools::install_local("GEOquery-main.zip")
#library(GEOquery)
#options(timeout=200000)
#getOption("timeout")
#options( 'download.file.method.GEOquery' = 'libcurl' )
#a = getGEO("GSE5859",destdir = ".")
#library(GSE5859)
#data(GSE5859)
#This is the original dataset from which we selected the subset used in GSE5859Subset.
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GEOquery")
options(stringsAsFactors = F)
# 注意查看下载文件的大小，检查数据 
f='GSE5859.Rdata'

library(GEOquery)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
options( 'download.file.method.GEOquery' = 'libcurl' )
options(stringsAsFactors = F)
# 注意查看下载文件的大小，检查数据 
f='GSE5859.Rdata'

library(GEOquery)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
options( 'download.file.method.GEOquery' = 'libcurl' )
options(timeout =200000)
if(!file.exists(f)){
  gset <- getGEO('GSE5859', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE5859.Rdata')  ## 载入数据
#We can extract the gene expression data and sample information table using the Bioconductor functions exprs and pData like this:
geneExpression = exprs(gset[[1]])
#8793*208
sampleInfo = pData(gset[[1]])
gset[[1]]@annotation

## get microarray platform, "GPL201", find corresponding microarray-SYMBOL transform library
## http://www.bio-info-trainee.com/1399.html, according to the website, 
## GPL201 corresponds to hgfocus 
library(hgfocus.db)

## this contains all id-transform capability
ls("package:hgfocus.db")
## dataframe, which has probeid -SYMBOL relationship
## Note that one gene may have multiple probes
ids=toTable(hgfocusSYMBOL)
#tmp = select(hgfocus.db, keys=featureNames(gset[[1]]), keytype="PROBEID",
#columns=c("CHR"，'CHRLOC','SYMBOL'))
#8379
SYMBOL = toTable(hgfocusSYMBOL)$symbol
PROBEID = toTable(hgfocusSYMBOL)$probe_id
chrtmp = select(hgfocus.db, keys=featureNames(gset[[1]]), keytype="PROBEID",
             columns=c("CHR"))
chr = chrtmp[which(chrtmp$PROBEID %in% PROBEID),]
CHRtmp = ifelse(is.na(chr$CHR),'NA',paste0('chr',chrtmp$CHR))
CHR = rep(0,length(cclc[,1]))
for( i in 1:length(cclc[,1])){CHR[i] = median(CHRtmp[which(chrtmp$PROBEID == cclc[i,1])])}
chrloctmp = select(hgfocus.db, keys=featureNames(gset[[1]]), keytype="PROBEID",
                columns=c("CHRLOC"))
#CHRLOC = toTable(hgfocusCHRLOC)$start_location
#chrloc = chrloctmp[match(chrloctmp$PROBEID,PROBEID),]
chrloc = chrloctmp[which(chrloctmp$PROBEID %in% PROBEID),]
CHRLOC = chrloc$CHRLOC #14872
cclc= count(chrloc,chrloc$PROBEID)
CHRLOC = rep(0,length(cclc[,1]))
for( i in 1:length(cclc[,1])){CHRLOC[i] = median(chrloc$CHRLOC[which(chrloc$PROBEID == cclc[i,1])])}
#geneAnnotation = data.frame(c(PROBEID,SYMBOL,CHR,CHRLOC),colnames = c('PROBEID','SYMBOL','CHR','CHRLOC'))
geneAnnotation = data.frame(c(PROBEID,SYMBOL,CHR,CHRLOC))
colnames(geneAnnotation) = c('PROBEID','SYMBOL','CHR','CHRLOC')
#8393
#PROBEID = toTable(hgfocusCHR)$probe_id
#chrtmp[which(chrtmp$probe_id == PROBEID)]
#chrtmp[match(chrtmp$probe_id,PROBEID)]
#CHR = ifelse(is.na(toTable(hgfocusCHR)$chromosome),'NA',paste0('chr',toTable(hgfocusCHR)$chromosome))
#
#CHRLOC = toTable(hgfocusCHRLOC)$start_location
#sampleInfo
GEOINFO = sampleInfo$geo_accession
updtdate = sampleInfo$last_update_date
#'''fdate<-function(x){
#  x = gsub(' ','-',x)
#  y = gsub('Jan','01',x)
#  y = gsub('Feb','02',x)
#  y = gsub('Mar','03',x)
#  y = gsub('Apr','04',x)
#  y = gsub('May','05',x)
#  y = gsub('Jun','06',x)
#  y = gsub('Jul','07',x)
#  y = gsub('Aug','08',x)
#  y = gsub('Sep','09',x)
#  y = gsub('Oct','10',x)
#  y = gsub('Nov','11',x)
#  y = gsub('Dec','12',x)
#  as.Date(c(y),format = '%m-%d-%Y')
#}'''
sapply(updtdate,fdate)
fdate<-function(x){
  as.Date(c(gsub('Dec','12',gsub('Nov','11',gsub('Oct','10',gsub('Sep','09',gsub('Aug','08',gsub('Jul','07',gsub('Jun','06',gsub('May','05',gsub('Apr','04',gsub('Mar','03',gsub('Feb','02',gsub('Jan','01',gsub(' ','-',x)))))))))))))),format = '%m-%d-%Y')
}
date = fdate(updtdate)

###!!!downloaded year only contains 25, 11, and no ethnicity thus change back to subset data

#查看expression 是否在20内，为log转换后的
head(geneExpression)
#Familiarize yourself with the sampleInfo table. Note that some samples were processed at different times. This is an extraneous variable and should not affect the values in geneExpression. However, as we have seen in previous analyses, it does appear to have an effect so we will explore this here.
load('GSE5859.rda')
geneExpression = exprs(e)
sampleInfo = pData(e)
#You can extract the year from each date like this:
year = format(sampleInfo$date,"%y")
#Note that ethnic group and year is almost perfectly confounded:
tby = table(year,sampleInfo$ethnicity)
#year ASN CEU HAN
#02   0  32   0
#03   0  54   0
#04   0  13   0
#05  80   3   0
#06   2   0  24
#For how many of these years do we have more than one ethnicity represented?
#1 '05'  
#Repeat the above exercise, but now, instead of year, consider the month as well. Specifically, instead of the year variable defined above use:
month.year = format(sampleInfo$date,"%m%y")
#For what proportion of these month.year values do we have more than one ethnicity represented?
table(month.year, sampleInfo$ethnicity)
#month.year ASN CEU HAN
#0103   0  27   0
#0203   0  18   0
#0205   0   2   0
#0303   0   7   0
#0403   0   2   0
#0404   0   4   0
#0406   0   0  24
#0505   2   0   0
#0605  22   1   0
#0606   2   0   0
#0705   6   0   0
#0805  12   0   0
#0905   6   0   0
#1002   0   2   0
#1005  16   0   0
#1102   0  15   0
#1104   0   7   0
#1105  13   0   0
#1202   0  15   0
#1204   0   2   0
#1205   3   0   0


#Perform a t-test (use rowttests) comparing CEU samples processed in 2002 to those processed in 2003. Then use the qvalue package to obtain q-values for each gene.
CEU0203 = match(which(year %in%  c('02','03')),which(sampleInfo$ethnicity == 'CEU'))
rttsceu0203 = rowttests(geneExpression[,CEU0203])
rttsceu0203$p.value
any(is.finite(na.omit(rttsceu0203$p.value)))
library(qvalue)
#How many genes have q-values < 0.05 ?
  
#What is the estimate of pi0 provided by qvalue:
  
#Now perform a t-test (use rowttests) comparing CEU samples processed in 2003 to those processed in 2004. Then use the qvalue package to obtain q-values for each gene. How many genes have q-values less than 0.05?
CEU03 = match(which(year == '03'),which(sampleInfo$ethnicity == 'CEU')) 
CEU04 = match(which(year == '04'),which(sampleInfo$ethnicity == 'CEU')) 
rttsceu0304 = t.tests(t(geneExpression[,CEU03]),t(geneExpression[,CEU04]))
rttsceu0304$p.value
library(qvalue)
qvalue(rttsceu0304)
#Now we are going to compare ethnicities as was done in the original publication in which these data were first presented. Use the qvalue function to compare the ASN population to the CEU population. Once again, use the qvalue function to obtain q-values.
#How many genes have q-values < 0.05?
qvalue(na.omit(abs(seq(rttsceu0304$conf.int[1],rttsceu0304$conf.int[2],by=0.005))),fdr.level=0.05, pi0.method="bootstrap", adj=1.2,lambda = 0.0095,pfdr=TRUE)$qvalues  

#!!!!!!confounding solving
#Over 80% of genes are called differentially expressed between ethnic groups. However, due to the confounding with processing date, we need to confirm these differences are actually due to ethnicity. This will not be easy due to the almost perfect confounding. However, above we noted that two groups were represented in 2005. Just like we stratified by majors to remove the “major effect” in our admissions example, here we can stratify by year and perform a t-test comparing ASN and CEU, but only for samples processed in 2005.
#How many genes have q-values < 0.05 ?
ASN05 = match(which(year == '05'),which(sampleInfo$ethnicity == 'ASN')) 
CEU05 = match(which(year == '05'),which(sampleInfo$ethnicity == 'CEU')) 
rttsasnceu05 = t.test(t(geneExpression[,ASN05]),t(geneExpression[,CEU05]))
qvalue(na.omit(abs(seq(rttsasnceu05$conf.int[1],rttsasnceu05$conf.int[2],by=0.005))),fdr.level=0.05, pi0.method="bootstrap", adj=1.2,lambda = 0.035,pfdr=TRUE)$qvalues  
#Notice the dramatic drop in the number of genes with q-value < 0.05 when we fix the year. However, the sample size is much smaller in this latest analysis which means we have less power
table(sampleInfo$ethnicity[c(ASN05,CEU05)])
#ASN CEU HAN 
#3  80   0 
#To provide a more balanced comparison, we repeat the analysis, but now taking 3 random CEU samples from 2002. Repeat the analysis above, but comparing the ASN from 2005 to three random CEU samples from 2002. Set the seed at 3, set.seed(3)
#How many genes have q-values < 0.05 ?
set.seed(3)
ASN05 = match(which(year == '05'),which(sampleInfo$ethnicity == 'ASN')) 
CEU02 = sample(match(which(year == '02'),which(sampleInfo$ethnicity == 'CEU')),3) 
rttsasnceu = t.test(t(geneExpression[,ASN05]),t(geneExpression[,CEU02]))
qvalue(na.omit(abs(seq(rttsasnceu$conf.int[1],rttsasnceu$conf.int[2],by=0.005))),fdr.level=0.05, pi0.method="bootstrap", adj=1.2,lambda = 0.025,pfdr=TRUE)$qvaluetable(sampleInfo$ethnicity[c(ASN05,CEU02)])
table(sampleInfo$ethnicity[c(ASN05,CEU02)])
