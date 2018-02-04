library(glmnet)
siggenes= read.csv("siggenes_2.csv",header = FALSE)
siggenes = as.matrix(siggenes)
sigtfs = read.csv("sigtfs_2.csv",header=FALSE)
sigtfs = as.matrix(sigtfs)

yexp =read.csv("2704_geneswith_peakResult_cos_AVERAGE_GPSmooth.csv",header=TRUE)
ydata = yexp[,8:22]
ym = as.matrix(ydata) 
ym.c = scale(ym,center=TRUE,scale=FALSE)
ymt  = t(ym)
ymt.c = scale(ymt,center=TRUE,scale=FALSE)
dim(ym.c)
dim(ymt.c)
message(mean(ym.c[,1]),"  ", var(ym.c[,1]))
message(mean(ymt.c[,1]),"  ", var(ymt.c[,1]))

pdata = read.csv("pmatrix_paper_raw.csv",header=FALSE)
pm = as.matrix(pdata)
dim(pm)
pm = scale(pm,center=TRUE,scale=TRUE)
message(mean(pm[,1]),"  ", var(pm[,1]))
dim(pm)


paperdata = read.csv("amatrix_rustadpaper.csv",header=FALSE)
rustad = as.matrix(paperdata)
am=rustad
dim(rustad)
library(sigmoid)
rm = t(rustad)
rm[abs(rm)<1]=0
rmpie = 1-tanh(rm)# -1 to 0 to 1

#write.csv(rm,"test.csv")

log = c()
k = 0
while (k<5){k = k + 1;

  an = matrix(0,132,1)
  x1 = t(pm)
  
  
  for (i in 1:2704){
    #print(i);
    
    y1 = as.matrix(ymt.c[,i]);
    
    
    # # use cross validation
    model1 = cv.glmnet(x1,y1,type.measure = "mse", alpha = 1.0, nfolds = 7);
    pred = coef(model1,s="lambda.min")
    
    # # use penalty as prior
    # pie = as.numeric(rmpie[,i]);
    # model1 = glmnet(x1,y1,penalty.factor = pie);
    # lmmin = min(model1$lambda)
    # pred = coef(model1,s=lmmin)
    
    
    pred = as.matrix(pred[-1,])
    an = cbind(an,pred);}
  
  an = an[,2:2705]
  an = t(an)
  adiff = max(abs(an- am))
  am = an

pn = matrix(0,132,1)
x2 = scale(am,center=TRUE,scale=TRUE)
for (j in 1:15){
  #print(j);
  y2 = as.matrix(ym.c[,j]);
  model2 = cv.glmnet(x2,y2,type.measure = "mse", alpha = 1.0, nfolds = 10);
  pred = coef(model2,s="lambda.min");
  pred = as.matrix(pred[-1,]);
  pn = cbind(pn,pred);}

pn = pn[,2:16]
pdiff = max(abs(pn-pm))
pm= pn


log[k] = paste( "cycle", k, "  adiff:" , adiff, "    pdiff:", pdiff)
print(log[k])
}
log

rownames(am) = siggenes
colnames(am) = sigtfs
rownames(rustad) = siggenes
colnames(rustad) = sigtfs

counts = apply(am,1,function (r) {sum(r!=0)})
hist(counts,breaks=30)
sigtfs[am["Rv3849",] != 0]  
siggenes[am[,"Rv0678--"] != 0]
am["Rv2101","Rv3260c-whiB2"]

counts.rustad = apply(rustad,1,function (r) {sum(abs(r)>=1)})
hist(counts.rustad,breaks=30)
#look for other TF interaction
sigtfs[abs(rustad["Rv3849",]) >= 1]
#look for other targets
siggenes[abs(rustad[,"Rv0678--"]) >= 1]
rustad["Rv2101","Rv3260c-whiB2"]


################# comparison rustad vs model


commonresult <- c()
for (i in 1:132){
tf = sigtfs[i]
modeltargets = siggenes[am[,tf] != 0]
modellength = length(modeltargets)
rustadtargets = siggenes[abs(rustad[,tf]) >= 1]
rustadlength = length(rustadtargets)
commonlength = length(intersect(modeltargets, rustadtargets))
commonpercent = (commonlength / rustadlength ) * 100
commonresult[i] = paste(tf, modellength, rustadlength, commonlength, round(commonpercent,2), "%")
}
commonresult


############ scatterplot of y vs yhat
write.csv(am,"model_interactionsmatrix_am.csv")


yhat = am %*% pm

library(ggplot2)
library(grid)
library(gridExtra)


p1 <- qplot(yhat[,1],ym.c[,1])
p2 <- qplot(yhat[,2],ym.c[,2])
p3 <- qplot(yhat[,3],ym.c[,3])
p4 <- qplot(yhat[,4],ym.c[,4])
p5 <- qplot(yhat[,5],ym.c[,5])
p6 <- qplot(yhat[,6],ym.c[,6])
p7 <- qplot(yhat[,7],ym.c[,7])
p8 <- qplot(yhat[,8],ym.c[,8])
p9 <- qplot(yhat[,9],ym.c[,9])
p10 <- qplot(yhat[,10],ym.c[,10])
p11 <- qplot(yhat[,11],ym.c[,11])
p12 <- qplot(yhat[,12],ym.c[,12])
p13 <- qplot(yhat[,13],ym.c[,13])
p14 <- qplot(yhat[,14],ym.c[,14])
p15 <- qplot(yhat[,15],ym.c[,15])
figure <- grid.arrange(p1,p2,p3,p4, p5,p6,p7,p8,p9,p10,p11,p12, p13,p14,p15, ncol=4, top="Optimization of TF-Target Interaction using LASSO-glmnet")



############ summary of TF interactions in amatrix
col = 0
colindx =c()
tfinteractions = c()
for (i in 1:132){
  interactions = sum(am[,i]!=0.0) 
  if ((interactions)==0) {col = col + 1 ; colindx = c(colindx,i)}
               tfinteractions =c(tfinteractions,interactions) }

message ("  ZEROS columns: amatrix: ", col )
sigtfs[colindx]
tfinteractions = as.numeric(tfinteractions)
tfsummary = data.frame(sigtfs,tfinteractions)
tfsumorder = tfsummary[order(-tfinteractions),]
tfsumorder


write.csv(tfsumorder,"tfsummary.csv")
write.csv(commonresult,"commonresult.csv")



