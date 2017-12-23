
dat <- read.csv("1_mean_expression_cos_all2948_sig.csv",header=FALSE,row.names=1)
head(dat)
newdata=dat

library(mclust)
mixclust <- Mclust(newdata,G=5:25)  #it will test between 5 to 25
summary(mixclust)
emclust_result <- newdata
clustresult <- as.data.frame(mixclust$classification)
colnames(clustresult) = 'cluster'
emclust_result$clust <- clustresult$cluster
head(emclust_result)
write.csv(emclust_result, "emclust_cos_mean_all_result.csv")

dat <- read.csv("emclust_cos_mean_all_result.csv",header=TRUE,row.names=1)
head(dat)
emclust_result=dat


show_graph1 = function(lfcs,cnum) {
  plot(t(lfcs[1,1:15]),type="l",ylim=c(-3,3),xlab="",ylab="",main=cnum)  
  for (i in 2:length(rownames(lfcs))) {
    lines(t(lfcs[i,1:15])) } 
  med = sapply(lfcs[,1:15],median,2)
  lines(med,col=2) }

par(mfrow=c(4,5))
show_graph1(emclust_result[emclust_result$clust==1,],1)
show_graph1(emclust_result[emclust_result$clust==2,],2)
show_graph1(emclust_result[emclust_result$clust==4,],4)
show_graph1(emclust_result[emclust_result$clust==3,],3)
show_graph1(emclust_result[emclust_result$clust==5,],5)
show_graph1(emclust_result[emclust_result$clust==6,],6)
show_graph1(emclust_result[emclust_result$clust==8,],8)
show_graph1(emclust_result[emclust_result$clust==9,],9)
show_graph1(emclust_result[emclust_result$clust==10,],10)
show_graph1(emclust_result[emclust_result$clust==11,],11)
show_graph1(emclust_result[emclust_result$clust==12,],12)
show_graph1(emclust_result[emclust_result$clust==13,],13)
show_graph1(emclust_result[emclust_result$clust==14,],14)
show_graph1(emclust_result[emclust_result$clust==15,],15)
show_graph1(emclust_result[emclust_result$clust==16,],16)
show_graph1(emclust_result[emclust_result$clust==17,],17)
show_graph1(emclust_result[emclust_result$clust==18,],18)
show_graph1(emclust_result[emclust_result$clust==19,],19)
newdata = emclust_result


disteuc = dist(dat,method="euclidean")
fiteuc = hclust(disteuc,method="ward.D2")
groupseuc = cutree(fiteuc,k=30)
newdata = as.data.frame(newdata)
newdata$clust = groupseuc


colnames(newdata) <- c("3hr" , "6.5hr" , "9hr" , "12hr" , "18.5hr" , "21hr" , "27hr" , "31hr" , "33hr" , "36hr" , "39.5hr" , "42hr" , "45.5hr" , "52hr" , "55hr","clust")



show_graph1 = function(lfcs,cnum) {
  clusttotal=length(rownames(lfcs))
  med = sapply(lfcs[,1:15],median,2)
  plot(med,type="l",ylim=c(-2,2),color ='red',xlab="",ylab="",main=paste("clust#",cnum,"/ ",clusttotal, "genes"))
  for (i in 2:length(rownames(lfcs))) {
    lines(t(lfcs[i,1:15])) } }

png(file="1_10_19clusters_meanexp.png",width=12.0,height=10.0,units='in', res=300)

show_graph1 = function(lfcs,cnum) {
clusttotal=length(rownames(lfcs))
med = sapply(lfcs[,1:15],mean,2)
std = sapply(lfcs[,1:15],sd,2)
plot(med,type="l",ylim=c(-3,3),col ='red',xlab="",ylab="",main=paste("#",cnum,"(",clusttotal,"genes )"))#,"/ ",clusttotal, "genes"))
for ( j in 1:dim(lfcs)[1]){
   r = lfcs[j,1:15]
  #if (r<= med +std && r>= med - std){
     lines(t(r))
   #}
}
lines(med ,col ='red',lwd=2)#,xlab="",ylab="",main=paste("clust#",cnum,"/ ",clusttotal, "genes"))
}


par(mfrow=c(4,5))
show_graph1(newdata[newdata$clust==1,],1)
show_graph1(newdata[newdata$clust==2,],2)
show_graph1(newdata[newdata$clust==3,],3)
show_graph1(newdata[newdata$clust==4,],4)
show_graph1(newdata[newdata$clust==5,],5)
show_graph1(newdata[newdata$clust==6,],6)
show_graph1(newdata[newdata$clust==7,],7)
show_graph1(newdata[newdata$clust==8,],8)
show_graph1(newdata[newdata$clust==9,],9)
show_graph1(newdata[newdata$clust==10,],10)
show_graph1(newdata[newdata$clust==11,],11)
show_graph1(newdata[newdata$clust==12,],12)
show_graph1(newdata[newdata$clust==13,],13)
show_graph1(newdata[newdata$clust==14,],14)
show_graph1(newdata[newdata$clust==15,],15)
show_graph1(newdata[newdata$clust==16,],16)
show_graph1(newdata[newdata$clust==17,],17)
show_graph1(newdata[newdata$clust==18,],18)
show_graph1(newdata[newdata$clust==19,],19)
dev.off()

png(file="6clusters_meanexp.png",width=12.0,height=7.0,units='in', res=300)
par(mfrow=c(2,3))
show_graph1(newdata[newdata$clust==2,],1)
show_graph1(newdata[newdata$clust==11,],2)
show_graph1(newdata[newdata$clust==1,],3)
show_graph1(newdata[newdata$clust==19,],4)
show_graph1(newdata[newdata$clust==10,],5)
show_graph1(newdata[newdata$clust==16,],6)
dev.off()
