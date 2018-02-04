
library('gplots')
library('ggplot2')
library('knitr')
library('reshape2')
library('RColorBrewer')
library('WGCNA')
library ('igraph')
library('flashClust')
library(lattice)
###FOR WGCNA
# source("http://bioconductor.org/biocLite.R") 
# biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
# install.packages("WGCNA")


raw_counts <- read.csv('1_mean_expression_cos_all2948_sig.csv', header=FALSE,row.names=1)
head(raw_counts)
raw_counts = t(raw_counts)


#######   correlation network using free scale topology


######    this combines correlation and eculidean distance
cordist <- function(dat) {
    cor_matrix  <- cor(t(dat))
    # calculate distance eculidean, show diagnola and upper triangle of matrix
    dist_matrix <- as.matrix(dist(dat, diag=TRUE, upper=TRUE))
    #log1p(x) function computes log(x+1) accurately.
    #dist_matrix <- log1p(dist_matrix)
    #bring between 0 and 1 
    dist_matrix <- 1 - (dist_matrix / max(dist_matrix))
    # sign for -1 or +1 - combine correlation and E(graph)
    sign(cor_matrix) * ((abs(cor_matrix) + dist_matrix)/ 2)
}



sim_matrix <- cordist(raw_counts)
dim(sim_matrix)
# transformation to push high/down correlation values to preseve strongest correlation
adj_matrix <- adjacency.fromSimilarity(sim_matrix, power=10, type='signed')

rm(sim_matrix)
gc()

graph <- graph_from_adjacency_matrix(adj_matrix, mode=c("undirected"), weighted = TRUE)

V(graph)
E(graph)
graph <- simplify(graph,remove.loops=T)

#plot(E(graph)$weight)
sum(E(graph)$weight >=.9)
sum(E(graph)$weight >=.8)
sum(E(graph)$weight >=.7)
sum(E(graph)$weight >=.6)
sum(E(graph)$weight >=.5)

#sum(degree(graph)<1)

#graph <- simplify(graph,remove.loops=T)
graph=delete.edges(graph, which(E(graph)$weight <.8)) # here's my condition.
graph=delete.vertices(graph,which(degree(graph)<1))


#write.graph(graph, "network_periodic.graphml", format='graphml')
write.graph(graph, "4_network_0.9.graphml", format='graphml')






########### correlation model using wcgna function

gene.names=colnames(raw_counts)
SubGeneNames=gene.names

powers = c(c(1:10), seq(from = 2, to=30, by=2))
sft=pickSoftThreshold(raw_counts,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "signed")

# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#####  returns sum of the adjacency to the other node
 ## get k
con = softConnectivity(raw_counts,corFnc = "cor", corOptions = "use ='p'",
type = "signed", power = 16,blockSize = 1500,minNSamples = NULL,verbose = 2, indent = 0)
## hist plot at k
hist(con)   
summary(con)
#plot(log(temp),log(temp^))

scaleFreePlot(con,nBreaks = 10,truncated = FALSE,removeFirst = FALSE,main = "")

t1 = scaleFreeFitIndex(con, nBreaks = 10, removeFirst = FALSE)

####

df= data.frame(x = seq(1:2947), y = sort(con))
head(df)
df$cutx <- as.numeric(cut(df$x, breaks = 10))
df.agg <- aggregate(y ~ cutx, data = df, mean)
plot(log10(df.agg$y),log10(df.agg$y ^ 16),type='l')
par(new=TRUE)
plot(log10(df.agg$y))
####


softPower = 30;

#calclute the adjacency matrix
adj= adjacency(raw_counts,type = "signed", power = softPower);

#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(raw_counts,networkType = "signed", TOMType = "signed", power = softPower);

colnames(TOM) =rownames(TOM) =SubGeneNames
dissTOM=1-TOM

#hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average");

#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="",cex=0.3);


# Set the minimum module size
minModuleSize = 25;

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")


#discard the unassigned genes, and focus on the rest
restGenes= (dynamicColors != "grey")
length(restGenes)
diss1=1-TOMsimilarityFromExpr(raw_counts[,restGenes], power = softPower)

#plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

colnames(diss1) =rownames(diss1) =SubGeneNames[restGenes]
hier1=flashClust(as.dist(diss1), method="average" )
#plotDendroAndColors(hier1, dynamicColors[restGenes], "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")


#set the diagonal of the dissimilarity to NA 
diag(diss1) = NA;

#Visualize the Tom plot. Raise the dissimilarity matrix to the power to bring out the module structure
sizeGrWindow(7,7)
TOMplot(diss1, hier1, as.character(dynamicColors[restGenes]))

## get modules from the gene

module_colors= setdiff(unique(dynamicColors), "grey")
modules=c()
for (color in module_colors){
  module=SubGeneNames[which(dynamicColors==color)]
  modules[color] = list(module)
  #write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}
length(modules)


### heat map of each module

dat = raw_counts[,modules$purple]

my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
col_breaks = c(seq(-0.99,-0.75,length=100),  # for red
               seq(-0.74,0.74,length=100),           # for yellow
               seq(0.75,1,length=100))
heatmap.2(as.matrix(cor(dat)),
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row" ,   # only draw a row dendrogram
          #Colv="none", Rowv="none",
          #cexRow=1, cexCol=1, lwid=c(0.1,1), lhei=c(0.1,1),
          #margins=c(5,25)
          key=TRUE) # column label for 10 and row label size 25

####################################3


genenames = modules[[11]]
plot(raw_counts[,genenames[1]],type="l",ylim=c(-2,2),xlab="",ylab="",main=length(genenames))  
for (j in 2:length(genenames)) {
  lines(raw_counts[,genenames[j]]) } 


par(mfrow=c(4,5))
for (i in 1:length(modules)){
  genenames = modules[[i]]
  print(length(genenames))
  plot(raw_counts[,genenames[1]],type="l",ylim=c(-2,2),xlab="",ylab="",main=length(genenames))  
  for (j in 2:length(genenames)) {
    lines(raw_counts[,genenames[j]]) } 
  med = apply(raw_counts[,genenames[]],1,median)
  lines(med,col=2)
  }

write.csv(t(raw_counts[,modules$purple]),'module.csv')

medians=c()
#par(mfrow=c(4,5))
for (i in 1:length(modules)){
  genenames = modules[[i]]
  med = apply(raw_counts[,genenames[]],1,median)
  #plot(med,col=2,type="l")
  medians[i] =  paste(med,collapse=",")
}
write.table(medians,"medians.csv",sep=",",quote = FALSE)




#### try straight correlation values

my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
col_breaks = c(seq(-0.99,-0.75,length=100),  # for red
               seq(-0.74,0.74,length=100),           # for yellow
               seq(0.75,1,length=100))
heatmap.2(as.matrix(cor(t(raw_counts))),
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row" ,   # only draw a row dendrogram
          #Colv="none", Rowv="none",
          #cexRow=1, cexCol=1, lwid=c(0.1,1), lhei=c(0.1,1),
          #margins=c(5,25)
          key=TRUE) # column label for 10 and row label size 25


#### try hard threshholding 

data = cor(t(raw_counts))
# temp = as.data.frame(apply(data, 2, function(x) sum(x)))
# temp = as.data.frame(apply(data, 2, function(x) if(abs(x)<0.9){x=0.0}))

data[abs(data) < .9] = 0

my_palette <- colorRampPalette(c("yellow",  "red"))(n = 199)
col_breaks = c(seq(0.0,0.8,length=100),  # for red
               seq(0.9,1.0,length=100))

heatmap.2(as.matrix(data),
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row" ,   # only draw a row dendrogram
          #Colv="none", Rowv="none",
          #cexRow=1, cexCol=1, lwid=c(0.1,1), lhei=c(0.1,1),
          #margins=c(5,25)
          key=TRUE) # column label for 10 and row label size 25





