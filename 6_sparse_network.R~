library(glmnet)
set.seed(1)

xexp =read.csv("1_sig_tf.csv",header=FALSE)
dim(xexp)
head(xexp)
xdata = xexp[,3:17]
head(xdata)
xm = as.matrix(xdata) 
xmt  = t(xm)
dim(xmt)
colnames(xmt) = paste(xexp[,1],'/',xexp[,2])
#xmt= scale(xmt,center=TRUE,scale=TRUE)
message(mean(xmt[,1]),"  ", var(xmt[,1]))



yexp =read.csv("1_mean_expression_cos_all4018_cellcycle.csv",header=FALSE)
head(yexp)
ydata = yexp[,3:17]
head(ydata)
ym = as.matrix(ydata) 
ymt  = t(ym)
dim(ymt)
colnames(ymt) = paste(yexp[,1],'/',yexp[,2])
#ymt= scale(ymt,center=TRUE,scale=TRUE)

result = c()

#xm = as.matrix(xmt[1:14,])
xm = as.matrix(xmt)
for (i in 1:dim(yexp)[1]){

  #y1 = as.matrix(ymt[-1,i])
   
  y1 = as.matrix(ymt[,i])
  
  #model1 = cv.glmnet(xm,y1,type.measure="mse",nfolds = 15 )#, alpha =1.0 ,standardize=FALSE)
  #lmin = model1$lambda.min
  #lmin = model1$lambda.1se ## this is too high
  
  
  model1 = glmnet(xm,y1)
  lmin = min(model1$lambda)
  
  pred = coef(model1,s=lmin)
  pred = as.matrix(pred[-1,])



  parents = which(pred!=0.0)
  coefs = pred[which(pred !=0)]

  # yhat = xm %*% pred
  # yhat_cpd = 0.0
  # for (j in 1:length(yhat)){
  #   #yhat_cpd = yhat_cpd * pnorm(y1[j],yhat[j],0.5)
  #   yhat_cpd = yhat_cpd + abs(abs(yhat[j])-abs(y1[j]))^2
  # }
  # cbind(y1,yhat,yhat_cpd)

  #result[i]=paste(i,",",length(parents),",", paste(parents,collapse=","), ",","coefs",",",paste(coefs,collapse=","), ",","cpd",",",paste(yhat_cpd,collapse=","))

  temp =c()
  for(k in 1:length(parents)){

    temp[k] = paste(parents[k],",",coefs[k])

  }

  result[i] = noquote(paste(i,",",length(parents),",", paste(temp,collapse=",")))
  #print(result[i])

  # result_prob[i] =paste(yhat_cpd)



}


# result_prob = as.numeric(result_prob)
# which(result_prob >0.0002)
# plot(result_prob)

write.table(result,"output_tf_parents_coefs_cellcycle.csv",row.names = FALSE, col.names = FALSE,sep=",",quote = FALSE)

#write.csv(result,"output_tf_parents_coefs_cellcycle2.csv",row.names = FALSE, col.names = FALSE,sep=",",quote = FALSE)


#write.csv(which(result_prob >0.0002),"output_tf_parents_probs_cellcycle_index.csv",row.names = FALSE)





#######check for best cut off for tfs
var=0
var_result =c()
while( var < 16){


  var_residuals = 0.0

for (i in 1:dim(yexp)[1]){

  y1 = as.matrix(ymt[-1,i])

  xm = as.matrix(xmt[1:14,-i])


  model1 = glmnet(xm,y1, alpha =1.0 ,standardize=FALSE,dfmax = var)

  lmmin = min(model1$lambda)
  pred = coef(model1,s=lmmin)
  pred = as.matrix(pred[-1,])



  parents = which(pred!=0.0)
  coefs = pred[which(pred !=0)]

  yhat = xm %*% pred
  yhat_cpd = 0.0
  for (j in 1:length(yhat)){
    #yhat_cpd = yhat_cpd * pnorm(y1[j],yhat[j],0.5)
    yhat_cpd = yhat_cpd + abs(abs(yhat[j])-abs(y1[j]))^2
  }
  #cbind(y1,yhat,yhat_cpd)

  #result[i]=paste(i,",",length(parents),",", paste(parents,collapse=","), ",","coefs",",",paste(coefs,collapse=","), ",","cpd",",",paste(yhat_cpd,collapse=","))

  # temp =c()
  # for(k in 1:length(parents)){
  #
  #   temp[k] = paste(parents[k],",",coefs[k])
  #
  # }

  # result[i] = paste(i,",",length(parents),",", paste(temp,collapse=","))
  #print(result[i])


  var_residuals= var_residuals + yhat_cpd


}

  var_result[var] = paste(var_residuals)

  var = var + 1
  }

var_result = as.numeric(var_result)

plot(var_result, ylab='Sum of squared residuals',xlab='Number of Transcription Factor in the model')


#####################################
var_result_diff = c() # apply(as.matrix(var_result), 1:length(var_result)-1, function(x) abs(var_result[x] - var_result[x+1]))
for (i in 1:length(var_result)){
  var_result_diff[i] = abs(var_result[i+1]- var_result[i])
  message(i)
}

plot(var_result_diff)



# x=1:15
# plot( x, ymt[,177], type="l", col="red" )
# par(new=TRUE)
# plot( x, ymt[,1], type="l", col="black" )
# par(new=TRUE)
# plot( x, ymt[,7], type="l", col="black" )
# par(new=TRUE)
# plot( x, ymt[,28], type="l", col="black" )
# par(new=TRUE)
# plot( x, ymt[,87], type="l", col="black" )



##################################

##################################
