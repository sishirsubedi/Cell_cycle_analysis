
yexp =read.csv("2704_geneswith_peakResult_cos_AVERAGE_GPSmooth.csv",header=TRUE)
ydata = yexp[,8:22]
ym = as.matrix(ydata) 
ym.c = scale(ym,center=TRUE,scale=TRUE)
ymt  = t(ym)
ymt.c = scale(ymt,center=TRUE,scale=TRUE)
dim(ym.c)
dim(ymt.c)
message(mean(ym.c[,1]),"  ", var(ym.c[,1]))
message(mean(ymt.c[,1]),"  ", var(ymt.c[,1]))

siggeneswn = as.matrix(paste(yexp[,1],"-",yexp[,2]))
siggenes = as.matrix(yexp[,1])
sigtfs = read.csv("sigtfs_132.csv",header=FALSE)
sigtfs = as.matrix(sigtfs)
sigtfswn = read.csv("sigtfs_wn.csv",header=FALSE)
sigtfswn = as.matrix(sigtfswn)


pdata = read.csv("pmatrix_paper_raw.csv",header=FALSE)
pm = as.matrix(pdata)
dim(pm)
pm.c = scale(pm,center=TRUE,scale=TRUE)
message(mean(pm.c[,1]),"  ", var(pm.c[,1]))


################## linear model

x1 = t(pm) 
x1 = scale(x1,center=TRUE,scale=TRUE)
y1 = as.matrix(ymt.c[,1]);
df = data.frame(y1,x1)
model1 = lm(df)
summary(model1)
# does not work becaue of many variables 
# because two or more of your independent variables
# are perfectly collinear 



## trying lasso method
library(glmnet)
x1 = t(pm)
y1 = as.matrix(ymt.c[,1]);
model1 = glmnet(x1,y1, alpha = 1.0)
summary(model1)
model1
lmmin = min(model1$lambda)
pred = coef(model1,s=lmmin)

## here pred gives sparse coefficients


################ PCA analysis on Pmatrix
pdata = read.csv("pmatrix_paper_raw.csv",header=FALSE)
pm = as.matrix(pdata)
pm = t(pm)
dim(pm)
pm = scale(pm,center=TRUE,scale=TRUE)
prin_comp <- prcomp(pm, scale. = F)
# mean to variables
prin_comp$center
# standard deviation of variables
prin_comp$scale
# rotation values
prin_comp$rotation
dim(prin_comp$rotation)
## rotated values
prin_comp$x
dim(prin_comp$x)
biplot(prin_comp, scale = 0)
std_dev <- prin_comp$sdev
pr_var <- std_dev^2
#proportion of variance explained
prop_varex <- pr_var/sum(pr_var)
plot(prop_varex, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")


pcx1 = prin_comp$x[,1:7]
y1 = as.matrix(ymt.c[,2000]);
df = data.frame(y1,pcx1)
model1 = lm(df)
summary(model1)

colnames(pm) = sigtfswn

###

t1 =  prin_comp$x %*% pm 
colnames(t1) = sigtfswn
write.csv(t1,"test1.csv")
###
an = matrix(0,7,1)
x1 = pcx1
for (i in 1:2704){
  #print(i);
  
  y1 = as.matrix(ymt.c[,i]);
  df = data.frame(y1,x1)
  model1 = lm(df)
  pred = as.matrix(coef(model1))
  pred = as.matrix(pred[-1,])
  pred = as.matrix(pred)
  an = cbind(an,pred);}

an = an[,2:2705]
am = t(an)


rownames(am) = siggenes

yhat = pcx1 %*%  an[,1:15] 
library(ggplot2)
library(grid)
library(gridExtra)


p1 <- qplot(yhat[,1],ymt.c[,1])
p2 <- qplot(yhat[,2],ymt.c[,2])
p3 <- qplot(yhat[,3],ymt.c[,3])
p4 <- qplot(yhat[,4],ymt.c[,4])
p5 <- qplot(yhat[,5],ymt.c[,5])
p6 <- qplot(yhat[,6],ymt.c[,6])
p7 <- qplot(yhat[,7],ymt.c[,7])
p8 <- qplot(yhat[,8],ymt.c[,8])
p9 <- qplot(yhat[,9],ymt.c[,9])
p10 <- qplot(yhat[,10],ymt.c[,10])
p11 <- qplot(yhat[,11],ymt.c[,11])
p12 <- qplot(yhat[,12],ymt.c[,12])
p13 <- qplot(yhat[,13],ymt.c[,13])
p14 <- qplot(yhat[,14],ymt.c[,14])
p15 <- qplot(yhat[,15],ymt.c[,15])
figure <- grid.arrange(p1,p2,p3,p4, p5,p6,p7,p8,p9,p10,p11,p12, p13,p14,p15, ncol=4, top="Observed vs PCA predicted for first 15 genes")


pcx1 = prin_comp$x[,1:7]
result = c()
for (i in 1:2704){
  gene = siggeneswn[i]
  y1 = as.matrix(ym[i,])
  corlist = vector()
  for (j in 1:7){
    y2 = as.matrix(pcx1[,j])
    corl = cor(y1,y2)
    corlist = c(corlist, corl)
  }
  result[i] = paste(gene,",",paste(corlist,collapse=","))
  }
write.csv(result,"pc_correlations_genes.csv")

pcx1 = prin_comp$x[,1:7]
result = c()
for (i in 1:132){
  gene = sigtfswn[i]
  y1 = as.matrix(pm[i,])
  corlist = vector()
  for (j in 1:7){
    y2 = as.matrix(pcx1[,j])
    corl = cor(y1,y2)
    corlist = c(corlist, corl)
  }
  result[i] = paste(gene,",",paste(corlist,collapse=","))
}
write.csv(result,"pc_correlations_tf.csv")

plot(pcx1[,1])




y1 = as.matrix(pm[1,])
  y2 = as.matrix(pcx1[,1])
  cor(y1,y2)
  cov(y1,y2)
  
  
  
# svd analysis
  pdata = read.csv("pmatrix_paper_raw.csv",header=FALSE)
  pm = as.matrix(pdata)
  pm = t(pm)
  dim(pm)
  pm = scale(pm,center=TRUE,scale=TRUE)
  svd1 = svd(pm)
  svd1
  plot(svd1$d)
  
