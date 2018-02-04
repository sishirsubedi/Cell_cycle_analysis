
library(lars)


yexp =read.csv("2704_geneswith_peakResult_cos_AVERAGE_GPSmooth.csv",header=TRUE)
ydata = yexp[,8:22]
ym = as.matrix(ydata) 
ym.c = ym # scale(ym,center=TRUE,scale=TRUE)
ymt  =t(ym)
ymt.c = ymt #scale(ymt,center=TRUE,scale=TRUE)
dim(ym.c)
dim(ymt.c)
message(mean(ym.c[,1]),"  ", var(ym.c[,1]))
message(mean(ymt.c[,1]),"  ", var(ymt.c[,1]))

adata = read.csv("amatrix_rustadpaper.csv",header=FALSE)
am = as.matrix(adata)
am = cbind(matrix(1,2704,1),am)
#am = scale(am,center=TRUE,scale=TRUE)
message(mean(am[,1]),"  ", var(am[,1]))
dim(am)



pdata = read.csv("pmatrix_paper_raw.csv",header=FALSE)
pm = as.matrix(pdata)
dim(pm)
pm = rbind(matrix(1,1,15),pm)
message(mean(pm[,1]),"  ", var(pm[,1]))
dim(pm)


k = 0
while(k<1){

	k = k + 1

	pn = matrix(0,133,1);
        x1=am
	for (i in 1:15){
		y1 = as.matrix(ym.c[,i]);
                model1 = lars(x1,y1,normalize=TRUE,intercept=FALSE);
                step1 = which.min(model1$Cp);
		larspd1 = predict(model1,x1,s=step1,type="coefficients",mode="step");
		coef1 = larspd1$coefficients;
 		coef1m = as.matrix(coef1);
    pn = cbind(pn,coef1m);
          }
    pn = pn[,2:16]
    pdiffs = pn -pm
    pm = pn

  message (" CYCLE: ",k, "  ZEROS: pmatrix: ", sum(pm == 0.0))
  message (" CYCLE: ",k, "  pdiffs: ", max(pdiffs))

  col = 0

  for (i in 1:15){
    if ((sum(pm[,i]==0.0))==133) {col = col + 1;}}
  message (" CYCLE: ",k, "  ZEROS columns: pmatrix: ", col )


	an = matrix(0,133,1)
	
	x2 = t(pm)
	
	
	for (j in 1:2704){
	  print(j)
		y2 = as.matrix(ymt.c[,j]);
    model2 = lars(x2,y2,normalize=TRUE,intercept=FALSE, use.Gram=TRUE);
    step2 = which.min(model2$Cp);
		larspd2 = predict(model2,x2,s=step2,type="coefficients",mode="step");
		coef2 = larspd2$coefficients;
 		coef2m = as.matrix(coef2);
    an = cbind(an,coef2m);}
  an = an[,2:2705]
  an2 = t(an)
	adiffs = an2 -am
  am = an2
  message (" CYCLE: ",k, "  ZEROS:  amatrix: " , sum(am == 0.0))
  message (" CYCLE: ",k, "  adiffs: " , max(adiffs))
  
  col = 0
  for (i in 1:133){
    if ((sum(am[,i]==0.0))==2704) {col = col + 1}}
  message (" CYCLE: ",k, "  ZEROS columns: amatrix: ", col )
 }
yhat = am %*% pm
yhat.c = scale(yhat,center=TRUE)
plot(ym.c[,1],yhat.c[,1])

y2 = as.matrix(ymt.c[,1]);
model2 = lars(x2,y2,normalize=FALSE,intercept=FALSE, use.Gram=TRUE);
step2 = which.min(model2$Cp);
larspd2 = predict(model2,x2,s=step2,type="coefficients",mode="step");
coef2 = larspd2$coefficients;
coef2m = as.matrix(coef2);

