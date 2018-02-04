

#problem with non fourier frequency 





# 
# library(MetaCycle)
# cycYeastCycle
# write.csv(cycYeastCycle, file="cycYeastCycle.csv", row.names=FALSE)
# cyc <- meta2d(infile="cycYeastCycle.csv",filestyle="csv",
#               minper=80, maxper=96, timepoints=seq(2, 162, by=16),
#               outputFile=FALSE, ARSdefaultPer=85, outRawData=TRUE )
# 
# write.csv(cyc, file="yeastcyc_output.csv", row.names=FALSE)
# 
# 
# cosMtb = read.csv("cos_smooth_divbymean_adjtp_eventp_2.csv", header = TRUE)
# coscyc <- meta2d(infile="cos_smooth_divbymean_adjtp_eventp_2.csv",filestyle="csv",
#                minper=24, maxper=30, timepoints=seq(3, 54, by=3),
#                outputFile=FALSE, ARSdefaultPer=27, outRawData=TRUE )
# 
#  write.csv(coscyc, file="cosMtb_output.csv", row.names=FALSE)
# 
# 
#  rvcyc <- meta2d(infile="rv_smooth_divbymean_adjtp_eventp_2.csv",filestyle="csv",
#                   minper=32, maxper=40, timepoints=seq(4, 60, by=4),
#                   outputFile=FALSE, ARSdefaultPer=36, outRawData=TRUE )
# 
#  write.csv(rvcyc, file="rvMtb_output.csv", row.names=FALSE)
# 

################## gene cycle


library(plyr)
# load GeneCycle library
library(GeneCycle)


cosraw1 <- read.csv("2_deseqnorm_cos1.csv",header=FALSE,row.names = 1)
cosraw2 <- read.csv("2_deseqnorm_cos2.csv",header=FALSE,row.names = 1)
rvraw1 <- read.csv("2_deseqnorm_rv1.csv",header=FALSE,row.names = 1)
rvraw2 <- read.csv("2_deseqnorm_rv2.csv",header=FALSE,row.names = 1)
dim(cosraw1)

cos = matrix(0,2948, 15)
rv = matrix(0,2948, 15)
for (x in 1:2948){
  for (y in 1:15){
    cos[x,y]= (cosraw1[x,y]+cosraw2[x,y])/2
    rv[x,y]= (rvraw1[x,y]+rvraw2[x,y])/2
  }
}
genes = rownames(cosraw1)
cosmat = as.matrix(cos)
cosmat[cosmat<0.05]=0.05
logcos = log10(cosmat)
cos = t(logcos)

rvmat = as.matrix(rv)
rvmat[rvmat<0.05]=0.05
logrv = log10(rvmat)
rv = t(logrv)

cosdf = t(dominant.freqs(cos, 1))
count(cosdf)
hist(cosdf)


rvdf = t(dominant.freqs(rv, 1))
count(rvdf)
hist(rvdf)


pval.cos <- fisher.g.test(cos)
fdr.cos <- fdrtool(pval.cos, statistic = "pvalue")
sum(fdr.cos$qval <0.05)
pv.cos = as.matrix(fdr.cos$qval)
rownames(pv.cos)=genes

pval.rv = fisher.g.test(rv)
fdr.rv = fdrtool(pval.rv, statistic = "pvalue")
sum(fdr.rv$qval <0.05)
pv.rv = as.matrix(fdr.rv$qval)
rownames(pv.rv)=genes

cosfinal = cbind(genes,pv.cos)
write.csv(cosfinal,"cosfinal.csv")
rvfinal = cbind(genes, pv.rv)
write.csv(rvfinal,"rvfinal.csv")


colnames(cos)= genes
pv.cos["Rv1295",1]
pv.cos["Rv1629",1]
plot(cos[1:15,"Rv1295"],type="l")
plot(cos[1:15,"Rv1629"],type="l")


t1 = avgp(cos[1:15,"Rv1295"],plot = TRUE)

########### TRY COS GP SMOOTH  ###########################
### problem is allmost all genes are periodic <0.05  #####


cosraw1 <- read.csv("cos_smooth_divbymean_adjtp_eventp_18.csv",header=TRUE,row.names = 1)
rvraw1 <- read.csv("rv_smooth_divbymean_adjtp_eventp_18.csv",header=TRUE,row.names = 1)
dim(cosraw1)


genes = rownames(cosraw1)

cosmat = as.matrix(cosraw1)
cosmat[cosmat<0.05]=0.05
logcos = log10(cosmat)
cos = t(logcos[,1:10])

rvmat = as.matrix(rvraw1)
rvmat[rvmat<0.05]=0.05
logrv = log10(rvmat)
rv = t(logrv)

cosdf = t(dominant.freqs(cos, 1))
count(cosdf)
hist(cosdf)


rvdf = t(dominant.freqs(rv, 1))
count(rvdf)
hist(rvdf)


pval.cos <- fisher.g.test(cos)
# fdr.cos <- fdrtool(pval.cos, statistic = "pvalue")
# sum(fdr.cos$qval <0.00001)
# pv.cos = as.matrix(fdr.cos$qval)
# rownames(pv.cos)=genes
sum(pval.cos<0.01)
pv.cos = as.matrix(pval.cos)
rownames(pv.cos)=genes



pval.rv = fisher.g.test(rv)
sum(pval.rv<0.01)
pv.rv = as.matrix(pval.rv)
rownames(pv.cos)=genes



colnames(cos)= genes

rvn= "Rv0222"
pv.cos[rvn,1]
plot(cos[,rvn],type="l")


rvn= "Rv1629"
pv.cos[rvn,1]
plot(cos[,rvn],type="l")


rvn= "Rv1295"
pv.cos[rvn,1]
plot(cos[,rvn],type="l")

rvn= "Rv1321"
pv.cos[rvn,1]
plot(cos[,rvn],type="l")

rvn= "Rv2181"
pv.cos[rvn,1]
plot(cos[,rvn],type="l")

rvn="Rv2181"
t = periodogram(cos[,rvn])
plot(t$spec)

cosfinal = cbind(genes,pv.cos)
write.csv(cosfinal,"cosfinal.csv")
rvfinal = cbind(genes, pv.rv)
write.csv(rvfinal,"rvfinal.csv")


############# end of GP Smooth test #################
#####################################################

library(plyr)
# load GeneCycle library
library("GeneCycle")
# load data set
data(caulobacter)
# how many samples and how many genes?
dim(caulobacter)
# average periodogram
avgp.caulobacter <- avgp(caulobacter, "Caulobacter")
avgp.caulobactor

avgp(caulobacter, "Caulobacter", plot=FALSE)




plot(cos[1:15,1254],type = "l")
# f1 = fft(cos[1:55,1254])
# magn = Mod(f1)
# magn.1 <- magn[1:(length(magn)/2)]
# plot(magn.1,type="l")
# x.axis <- 1:length(magn.1)/time
# plot(x=x.axis,y=magn.1,type="l")

# f1inv = fft(f1, inverse = TRUE)
# magninv = Mod(f1inv)

# plot(magninv,type="l", ylim = c(-3,20))
# par(new=TRUE)
# plot(cos[1:55,1254],type = "l", ylim = c(-3,100)) 

plot(cos[1:30,1200])
avgp.cos = avgp(cos[1:30,1254],"cos")
avgp.cos
maxindex = which.max(avgp.cos$avg.spec)
maxfreq = avgp.cos$freq[maxindex]
maxfreq
avgp.cos$freq




result = vector()
r = vector()

for ( i in 1:2948){
avgp.cos = avgp(cos[1:15,i],"cos",plot=FALSE)
maxindex = which.max(avgp.cos$avg.spec)
maxfreq = avgp.cos$freq[maxindex]
pv = fisher.g.test(cos[1:15,i])
result[i]= paste(genes[i],",",maxfreq,",",pv)
r = append(r,maxfreq)
}

count(r)
hist(r)


pval.cos <- fisher.g.test(cos)
fdr.cos <- fdrtool(pval.cos, statistic = "pvalue")
sum(fdr.cos$lfdr <0.01)

pval.rv = fisher.g.test(rv)
fdr.rv = fdrtool(pval.rv, statistic = "pvalue")
sum(fdr.rv$lfdr <0.01)

rv="Rv1629"
t = periodogram(cos[1:18,rv])
plot(t$spec)



r[1];r[100];r[265]
plot(cos[1:18,1])  # 0.11
plot(cos[1:18,100])  # 0.05
plot(cos[1:18,1254])      #0.16

avgp.cos = avgp(cos[1:15,2939],"cos")


write.csv(result,"result.csv")





set.seed(101)
acq.freq <- 200
time     <- 1
w        <- 2*pi/time
ts       <- seq(0,time,1/acq.freq)
trajectory <- 3*rnorm(101) + 3*sin(3*w*ts) + 2*(sin(10*w*ts))+10
plot(trajectory, type="l")
avgp.traj= avgp(trajectory,"cos")
which.max(avgp.traj$avg.spec)




plot(cos[1:15,1950],type = "l")
avgp.cos = avgp(cos[1:15,1254],"cos")
avgp.cos
maxindex = which.max(avgp.cos$avg.spec)
maxfreq = avgp.cos$freq[maxindex]
maxfreq


a = rnorm(n=100)
x = seq(1:100)
b = sin(2*2*3.14*x/100)
c = sin(1.5*2*3.14*x/100)
af = fft(a)
bf = fft(b)
cf=fft(c)
plot(Mod(af))
plot(Mod(bf))
plot(Mod(cf))
plot(Re(af)^2)
plot(Re(bf)^2)
plot(Re(cf)^2)


## Slow Discrete Fourier Transform (DFT) - e.g., for checking the formula
dft1 <- function(z) {
  n <- length(z)
  if(n == 0) return(z)
  k <- 0:(n-1)
  angle <- -1 * 2*pi * 1i * k/n
  vapply(1:n, function(h) sum(z * exp(angle*(h-1))), complex(1))
}

dft2 <- function(z) {
  n <- length(z)
  if(n == 0) return(z)
  k <- 0:(n-1)
  angle <- 1 * 2*pi * k/n
  vapply(1:n, function(h) sum(z * (cos(angle*(h-1))- (1i * sin(angle*(h-1))))), complex(1))
}


dfti <- function(z) {
  n <- length(z)
  if(n == 0) return(z)
  k <- 0:(n-1)
  angle <- -1 * 2*pi * 1i * k/n
  vapply(1:n, function(h) sum(z * exp(angle*(h-1))), complex(1))
}
  
dft2(c)
dfti(dft2(c))

####################
library(ptest)
cosraw1 <- read.csv("cos_smooth_divbymean_adjtp_eventp_18.csv",header=TRUE,row.names = 1)
rvraw1 <- read.csv("rv_smooth_divbymean_adjtp_eventp_18.csv",header=TRUE,row.names = 1)
dim(cosraw1)


genes = rownames(cosraw1)

result = c()



rvname = genes[1]
rvdata = cosraw1[rvname,]
fit = lm(t(rvdata)~x)
res = fit$residuals
plot(res)

x = seq(1:18)
for ( i in 1:2948){
  rvname = genes[i]
  rvdata = cosraw1[rvname,]
  fit = lm(t(rvdata)~x)
  res = fit$residuals
  result[i] = paste(res,collapse = ",")
}
write.csv(result,"residuals.csv")

cosraw1 <- read.csv("residuals.csv",header=FALSE,row.names = 1)
genes = rownames(cosraw1)
cosmat = as.matrix(cosraw1)
cos= t(cosmat)
colnames(cos) = genes


cosmat[cosmat<0.05]=0.05
logcos = log10(cosmat)
cos = t(logcos)
colnames(cos) = genes

rvmat = as.matrix(rvraw1)
rvmat[rvmat<0.05]=0.05
logrv = log10(rvmat)
rv = t(logrv)
colnames(rv)= genes



rv = cos[,"Rv0554"]
plot(rv,type="l")
t1 = ptestReg(rv,method="LS")
t1
t2 = ptestg(rv,method="Fisher")
t2
######################################## extended fisher test  ############
library(ptest)
cosraw1 <- read.csv("cos_smooth_divbymean_adjtp_eventp_18.csv",header=TRUE,row.names = 1)
dim(cosraw1)
genes = as.data.frame(rownames(cosraw1))
cosmat = as.matrix(cosraw1)

periodic <- read.table("periodic_398.txt",header=FALSE)

cosmat[cosmat<0.05]=0.05
logcos = log10(cosmat)
cos = t(logcos)
colnames(cos) = genes

result = vector()
r2 = vector()

for ( i in 1:2948){
  #print(i)
  rvname = genes[i,1]
  rv =  cos[,rvname]
  t1 = ptestg(rv,method="extendedRobust")
  maxfreq = t1$freq
  pv = t1$pvalue
  result[i]= paste(rvname,",",maxfreq,",",pv)
  r2 = append(r2,pv)
}


write.csv(result,"result_cos_residual_extendedRobust.csv")


###################### do all test to check if any method gives 398 genes

library(ptest)
cosraw1 <- read.csv("cos_smooth_divbymean_adjtp_eventp_18.csv",header=TRUE,row.names = 1)
dim(cosraw1)
genes = as.data.frame(rownames(cosraw1))
cosmat = as.matrix(cosraw1)

periodic <- read.table("periodic_398.txt",header=FALSE)

cosmat[cosmat<0.05]=0.05
logcos = log10(cosmat)
cos = t(logcos)
colnames(cos) = genes

result = vector()
r = vector()

ct1 = 0;ct2=0;ct3=0
ct4 = 0;ct5=0;ct6=0
for ( i in 1:2948){
  #print(i)
  rvname = genes[i,1]
  rv =  cos[,rvname]
  
  t1 = ptestg(rv,method="Fisher")
  pv1 = t1$pvalue
  if (pv1<0.01 &&  is.element(rvname,periodic[,1]) ){
   ct1 = ct1 + 1
  }

  
  t2 = ptestg(rv,method="robust")
  pv2 = t2$pvalue
  if (pv2<0.01 &&  is.element(rvname,periodic[,1]) ){
    ct2 = ct2 + 1
  }
  
  
  
  t3 = ptestg(rv,method="extended")
  pv3 = t3$pvalue
  if (pv3<0.01 &&  is.element(rvname,periodic[,1]) ){
    ct3 = ct3 + 1
  }
  
  
  t4 = ptestg(rv,method="extendedRobust")
  pv4 = t4$pvalue
  if (pv4<0.01 &&  is.element(rvname,periodic[,1]) ){
    ct4 = ct4 + 1
  }
  
  t5 = ptestReg(rv,method="LS")
  pv5 = t5$pvalue
  if (pv5<0.01 &&  is.element(rvname,periodic[,1]) ){
    ct5 = ct5 + 1
 
  }
  
  t6 = ptestReg(rv,method="LS")
  pv6 = t6$pvalue
  if (pv6<0.01 &&  is.element(rvname,periodic[,1]) ){
    ct6 = ct6 + 1
    
  }
  
  print(paste(i,ct1,ct2,ct3,ct4,ct5,ct6))

}



############

rvraw1 <- read.csv("rv_smooth_divbymean_adjtp_eventp_18.csv",header=TRUE,row.names = 1)
genes = rownames(rvraw1)
rvmat = as.matrix(rvraw1)
rvmat[rvmat<0.05]=0.05
logrv = log10(rvmat)
rv = t(logrv)
colnames(rv)= genes
result = vector()
r = vector()

for ( i in 1:2948){
  #print(i)
  rvname = genes[i]
  rvn =  rv[,rvname]
  t1 = ptestReg(rvn,method="LS")
  maxfreq = t1$freq
  pv = t1$pvalue
  result[i]= paste(rvname,",",maxfreq,",",pv)
  r = append(r,maxfreq)
}

write.csv(result,"result_rv.csv")


