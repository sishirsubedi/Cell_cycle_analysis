
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq")


library("DESeq")

library(ggplot2)
library(grid)
library(gridExtra)
library(reshape)

library(lattice) # correlation matrix

library(gplots) # for heatmap.2


############### READ sync files ##################

sync1 <- read.table("sync1.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_1", "rpkm_1"))
sync2 <- read.table("sync2.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_2", "rpkm_2"))
sync3 <- read.table("sync3.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_3", "rpkm_3"))
sync4 <- read.table("sync4.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_4", "rpkm_4"))
sync5 <- read.table("sync5.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_5", "rpkm_5"))
sync6 <- read.table("sync6.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_6", "rpkm_6"))
sync7 <- read.table("sync7.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_7", "rpkm_7"))
sync8 <- read.table("sync8.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_8", "rpkm_8"))
sync9 <- read.table("sync9.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_9", "rpkm_9"))
sync10 <- read.table("sync10.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_10", "rpkm_10"))
sync11 <- read.table("sync11.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_11", "rpkm_11"))
sync12 <- read.table("sync12.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_12", "rpkm_12"))
sync13 <- read.table("sync13.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_13", "rpkm_13"))
sync14 <- read.table("sync14.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_14", "rpkm_14"))
sync15 <- read.table("sync15.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_15", "rpkm_15"))
sync16 <- read.table("sync16.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_16", "rpkm_16"))


sync17 <- read.table("sync17.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_17","rpkm_17" ))
sync18 <- read.table("sync18.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_18","rpkm_18" ))
sync19 <- read.table("sync19.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_19","rpkm_19" ))
sync20 <- read.table("sync20.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_20","rpkm_20" ))
sync21 <- read.table("sync21.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_21","rpkm_21" ))
sync22 <- read.table("sync22.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_22","rpkm_22" ))
sync23 <- read.table("sync23.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_23","rpkm_23" ))
sync24 <- read.table("sync24.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_24","rpkm_24" ))
sync25 <- read.table("sync25.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_25","rpkm_25" ))
sync26 <- read.table("sync26.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_26","rpkm_26" ))
sync27 <- read.table("sync27.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_27","rpkm_27" ))
sync28 <- read.table("sync28.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_28","rpkm_28" ))
sync29 <- read.table("sync29.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_29","rpkm_29" ))
sync30 <- read.table("sync30.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_30","rpkm_30" ))
sync31 <- read.table("sync31.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_31","rpkm_31" ))
sync32 <- read.table("sync32.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_32","rpkm_32" ))




sync33 <- read.table("sync33.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_33", "rpkm_33"))
sync34 <- read.table("sync34.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_34", "rpkm_34"))
sync35 <- read.table("sync35.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_35", "rpkm_35"))
sync36 <- read.table("sync36.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_36", "rpkm_36"))
sync37 <- read.table("sync37.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_37", "rpkm_37"))
sync38 <- read.table("sync38.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_38", "rpkm_38"))
sync39 <- read.table("sync39.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_39", "rpkm_39"))
sync40 <- read.table("sync40.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_40", "rpkm_40"))
sync41 <- read.table("sync41.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_41", "rpkm_41"))
sync42 <- read.table("sync42.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_42", "rpkm_42"))
sync43 <- read.table("sync43.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_43", "rpkm_43"))
sync44 <- read.table("sync44.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_44", "rpkm_44"))
sync45 <- read.table("sync45.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_45", "rpkm_45"))
sync46 <- read.table("sync46.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_46", "rpkm_46"))
sync47 <- read.table("sync47.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_47", "rpkm_47"))
sync48 <- read.table("sync48.expr", header=FALSE, row.names = 1, col.names= c("rv","gname", "length","strand","count_48", "rpkm_48"))



sync49 <- read.table("sync49.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_49","rpkm_49" ))
sync50 <- read.table("sync50.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_50","rpkm_50" ))
sync51 <- read.table("sync51.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_51","rpkm_51" ))
sync52 <- read.table("sync52.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_52","rpkm_52" ))
sync53 <- read.table("sync53.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_53","rpkm_53" ))
sync54 <- read.table("sync54.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_54","rpkm_54" ))
sync55 <- read.table("sync55.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_55","rpkm_55" ))
sync56 <- read.table("sync56.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_56","rpkm_56" ))
sync57 <- read.table("sync57.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_57","rpkm_57" ))
sync58 <- read.table("sync58.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_58","rpkm_58" ))
sync59 <- read.table("sync59.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_59","rpkm_59" ))
sync60 <- read.table("sync60.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_60","rpkm_60" ))
sync61 <- read.table("sync61.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_61","rpkm_61" ))
sync62 <- read.table("sync62.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_62","rpkm_62" ))
sync63 <- read.table("sync63.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_63","rpkm_63" ))
sync64 <- read.table("sync64.expr", row.names =1, col.names = c("rvnum","gene", "length","orient","count_64","rpkm_64" ))

############################# normalize using deseq on count data - skip 0 hr #########

totalcount <- data.frame(sync2[4],sync3[4],sync4[4],sync5[4],sync6[4],sync7[4],sync8[4],sync9[4],sync10[4],sync11[4],sync12[4],sync13[4],sync14[4],sync15[4],sync16[4],sync18[4],sync19[4],sync20[4],sync21[4],sync22[4],sync23[4],sync24[4],sync25[4],sync26[4],sync27[4],sync28[4],sync29[4],sync30[4],sync31[4],sync32[4],sync34[4],sync35[4],sync36[4],sync37[4],sync38[4],sync39[4],sync40[4],sync41[4],sync42[4],sync43[4],sync44[4],sync45[4],sync46[4],sync47[4],sync48[4],sync50[4],sync51[4],sync52[4],sync53[4],sync54[4],sync55[4],sync56[4],sync57[4],sync58[4],sync59[4],sync60[4],sync61[4],sync62[4],sync63[4],sync64[4])
head(totalcount)

condition = factor (c("3hr" , "6.5hr" , "9hr" , "12hr" , "18.5hr" , "21hr" , "27hr" , "31hr" , "33hr" , "36hr" , "39.5hr" , "42hr" , "45.5hr" , "52hr" , "55hr","3hr" , "6.5hr" , "9hr" , "12hr" , "18.5hr" , "21hr" , "27hr" , "31hr" , "33hr" , "36hr" , "39.5hr" , "42hr" , "45.5hr" , "52hr" , "55hr", "3hr" , "6.5hr" , "9hr" , "12hr" , "18.5hr" , "21hr" , "27hr" , "31hr" , "33hr" , "36hr" , "39.5hr" , "42hr" , "45.5hr" , "52hr" , "55hr", "3hr" , "6.5hr" , "9hr" , "12hr" , "18.5hr" , "21hr" , "27hr" , "31hr" , "33hr" , "36hr" , "39.5hr" , "42hr" , "45.5hr" , "52hr" , "55hr"))


sample <- totalcount
sample[1:60] <- sapply(sample[1:60], as.integer)
cds = newCountDataSet( sample, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
totaldeseq <- counts( cds, normalized=TRUE )

colnames(totaldeseq) <- c( "deseq_2" ,"deseq_3", "deseq_4" , "deseq_5" ,"deseq_6" , "deseq_7" , "deseq_8" , "deseq_9" , "deseq_10" , "deseq_11" , "deseq_12" , "deseq_13" , "deseq_14" , "deseq_15" , "deseq_16", "deseq_18" ,"deseq_19", "deseq_20" , "deseq_21" ,"deseq_22" , "deseq_23" , "deseq_24" , "deseq_25" , "deseq_26" , "deseq_27" , "deseq_28" , "deseq_29" , "deseq_30" , "deseq_31" , "deseq_32",
"deseq_34" ,"deseq_35", "deseq_36" , "deseq_37" ,"deseq_38" , "deseq_39" , "deseq_40" , "deseq_41" , "deseq_42" , "deseq_43" , "deseq_44" , "deseq_45" , "deseq_46" , "deseq_47" , "deseq_48" ,"deseq_50" ,"deseq_51", "deseq_52" , "deseq_53" ,"deseq_54" , "deseq_55" , "deseq_56" , "deseq_57" , "deseq_58" , "deseq_59" , "deseq_60" , "deseq_61" , "deseq_62" , "deseq_63" , "deseq_64" )

head(totaldeseq )


#######################################################################################################

# look at separate rpkm for global gene expression profile

cos_1_rpkm <- data.frame(sync1[5],sync2[5],sync3[5],sync4[5],sync5[5],sync6[5],sync7[5],sync8[5],sync9[5],sync10[5],sync11[5],sync12[5],sync13[5],sync14[5],sync15[5],sync16[5])

cos_2_rpkm  <- data.frame(sync17[5],sync18[5],sync19[5],sync20[5],sync21[5],sync22[5],sync23[5],sync24[5],sync25[5],sync26[5],sync27[5],sync28[5],sync29[5],sync30[5],sync31[5],sync32[5])

rv_1_rpkm  <- data.frame(sync33[5],sync34[5],sync35[5],sync36[5],sync37[5],sync38[5],sync39[5],sync40[5],sync41[5],sync42[5],sync43[5],sync44[5],sync45[5],sync46[5],sync47[5],sync48[5])

rv_2_rpkm  <- data.frame(sync49[5],sync50[5],sync51[5],sync52[5],sync53[5],sync54[5],sync55[5],sync56[5],sync57[5],sync58[5],sync59[5],sync60[5],sync61[5],sync62[5],sync63[5],sync64[5])





cos_1_count <- data.frame(sync1[4],sync2[4],sync3[4],sync4[4],sync5[4],sync6[4],sync7[4],sync8[4],sync9[4],sync10[4],sync11[4],sync12[4],sync13[4],sync14[4],sync15[4],sync16[4])

cos_2_count <- data.frame(sync17[4],sync18[4],sync19[4],sync20[4],sync21[4],sync22[4],sync23[4],sync24[4],sync25[4],sync26[4],sync27[4],sync28[4],sync29[4],sync30[4],sync31[4],sync32[4])

rv_1_count <- data.frame(sync33[4],sync34[4],sync35[4],sync36[4],sync37[4],sync38[4],sync39[4],sync40[4],sync41[4],sync42[4],sync43[4],sync44[4],sync45[4],sync46[4],sync47[4],sync48[4])

rv_2_count <- data.frame(sync49[4],sync50[4],sync51[4],sync52[4],sync53[4],sync54[4],sync55[4],sync56[4],sync57[4],sync58[4],sync59[4],sync60[4],sync61[4],sync62[4],sync63[4],sync64[4])



cos_1_rpkm_sortbysync16 <- cos_1_rpkm[order(cos_1_rpkm$rpkm_16),]
tail(cos_1_rpkm_sortbysync16)



####### cos 1

cos_1 <- cos_1_count

p1<- qplot(rownames(cos_1),cos_1$count_1, data = cos_1, xlab="cos1", ylab="count") + theme(axis.text.x = element_text(size = 0)) 
p2<- qplot(rownames(cos_1),cos_1$count_2, data = cos_1, xlab="cos2", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p3<- qplot(rownames(cos_1),cos_1$count_3, data = cos_1, xlab="cos3", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p4<- qplot(rownames(cos_1),cos_1$count_4, data = cos_1, xlab="cos4", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p5<- qplot(rownames(cos_1),cos_1$count_5, data = cos_1, xlab="cos5", ylab="count") + theme(axis.text.x = element_text(size = 0))
p6<- qplot(rownames(cos_1),cos_1$count_6, data = cos_1, xlab="cos6", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p7<- qplot(rownames(cos_1),cos_1$count_7, data = cos_1, xlab="cos7", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p8<- qplot(rownames(cos_1),cos_1$count_8, data = cos_1, xlab="cos8", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p9<- qplot(rownames(cos_1),cos_1$count_9, data = cos_1, xlab="cos9", ylab="count") + theme(axis.text.x = element_text(size = 0))
p10<- qplot(rownames(cos_1),cos_1$count_10, data = cos_1, xlab="cos10", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p11<- qplot(rownames(cos_1),cos_1$count_11, data = cos_1, xlab="cos11", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p12<- qplot(rownames(cos_1),cos_1$count_12, data = cos_1, xlab="cos12", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p13<- qplot(rownames(cos_1),cos_1$count_13, data = cos_1, xlab="cos13", ylab="count") + theme(axis.text.x = element_text(size = 0))
p14<- qplot(rownames(cos_1),cos_1$count_14, data = cos_1, xlab="cos14", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p15<- qplot(rownames(cos_1),cos_1$count_15, data = cos_1, xlab="cos15", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p16<- qplot(rownames(cos_1),cos_1$count_16, data = cos_1, xlab="cos16", ylab=" ") + theme(axis.text.x = element_text(size = 0))


figure<- grid.arrange(p1,p2,p3,p4, p5,p6,p7,p8,p9,p10,p11,p12, p13,p14,p15,p16, ncol=4, top="")
ggsave("cos_1_global_exp_profile.png",figure,dpi=300)


#####################

cos_1 <- cos_1_rpkm

p1<- qplot(rownames(cos_1),cos_1$rpkm_1, data = cos_1, xlab="cos1", ylab="rpkm") + theme(axis.text.x = element_text(size = 0)) 
p2<- qplot(rownames(cos_1),cos_1$rpkm_2, data = cos_1, xlab="cos2", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p3<- qplot(rownames(cos_1),cos_1$rpkm_3, data = cos_1, xlab="cos3", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p4<- qplot(rownames(cos_1),cos_1$rpkm_4, data = cos_1, xlab="cos4", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p5<- qplot(rownames(cos_1),cos_1$rpkm_5, data = cos_1, xlab="cos5", ylab="rpkm") + theme(axis.text.x = element_text(size = 0))
p6<- qplot(rownames(cos_1),cos_1$rpkm_6, data = cos_1, xlab="cos6", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p7<- qplot(rownames(cos_1),cos_1$rpkm_7, data = cos_1, xlab="cos7", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p8<- qplot(rownames(cos_1),cos_1$rpkm_8, data = cos_1, xlab="cos8", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p9<- qplot(rownames(cos_1),cos_1$rpkm_9, data = cos_1, xlab="cos9", ylab="rpkm") + theme(axis.text.x = element_text(size = 0))
p10<- qplot(rownames(cos_1),cos_1$rpkm_10, data = cos_1, xlab="cos10", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p11<- qplot(rownames(cos_1),cos_1$rpkm_11, data = cos_1, xlab="cos11", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p12<- qplot(rownames(cos_1),cos_1$rpkm_12, data = cos_1, xlab="cos12", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p13<- qplot(rownames(cos_1),cos_1$rpkm_13, data = cos_1, xlab="cos13", ylab="rpkm") + theme(axis.text.x = element_text(size = 0))
p14<- qplot(rownames(cos_1),cos_1$rpkm_14, data = cos_1, xlab="cos14", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p15<- qplot(rownames(cos_1),cos_1$rpkm_15, data = cos_1, xlab="cos15", ylab=" ") + theme(axis.text.x = element_text(size = 0))
p16<- qplot(rownames(cos_1),cos_1$rpkm_16, data = cos_1, xlab="cos16", ylab=" ") + theme(axis.text.x = element_text(size = 0))


figure<- grid.arrange(p1,p2,p3,p4, p5,p6,p7,p8,p9,p10,p11,p12, p13,p14,p15,p16, ncol=4, top="", bottom="X-axis: Gene[Rv0001-Rv3924c]")
ggsave("1_figure_cos_1_rpkm_global_exp_profile.png",figure)

################################## plot last four

cos_1_rpkm_sortbysync16 <- cos_1_rpkm[order(cos_1_rpkm$rpkm_16),]
tail(cos_1_rpkm_sortbysync16)

##############################################################################3
##whisker plot



#totalcount <- as.data.frame(totalcount)
#totalrpkm <- as.data.frame(totalrpkm)
totaldeseq <- as.data.frame(totaldeseq)
#totaltmm <- as.data.frame(totaltmm)


totaldeseq_melt <- melt(totaldeseq)
totaldeseq_melt[,2] <- log2(totaldeseq_melt[,2] + 1)
figure <- ggplot(totaldeseq_melt, aes(x = variable, y = value, fill=variable)) + geom_boxplot() + xlab("") +ylab(expression(log[2](count + 1))) +ggtitle("") + theme(legend.position="none")
ggsave("3_figure_whiskerplot_deseqnormalization.png",figure,dpi=300)


############################################

PC=10
lfcdeseq = log((totaldeseq[,1:15]+totaldeseq[,16:30]+2*PC)/(totaldeseq[,31:45]+totaldeseq[,46:60]+2*PC),2)

colnames(lfcdeseq) <- c( "cos2" ,"cos3", "cos4" , "cos5" ,"cos6" , "cos7" , "cos8" , "cos9" , "cos10" , "cos11" , "cos12" , "cos13" , "cos14" , "cos15" , "cos16")
                         




df <- melt(lfcdeseq)
figure <- ggplot(df, aes(value, fill = variable)) + geom_density(alpha = 0.2)+ggtitle("")+
labs(x = "Log2 Fold Change")
ggsave("4_figure_densityplot_lfc.png",figure,dpi=300)


lfcrpkm = log((totalrpkm[,1:15]+totalrpkm[,16:30]+2*PC)/(totalrpkm[,31:45]+totalrpkm[,46:60]+2*PC),2)

lfccount = log((totalcount[,1:15]+totalcount[,16:30]+2*PC)/(totalcount[,31:45]+totalcount[,46:60]+2*PC),2)
###############################################################################


######### correlation matrix cos1 and Rv 1


cos_1_count <- data.frame(sync1[4],sync2[4],sync3[4],sync4[4],sync5[4],sync6[4],sync7[4],sync8[4],sync9[4],sync10[4],sync11[4],sync12[4],sync13[4],sync14[4],sync15[4],sync16[4])
colnames(cos_1_count) <- c( "deseq_1" ,"deseq_2" ,"deseq_3", "deseq_4" , "deseq_5" ,"deseq_6" , "deseq_7" , "deseq_8" , "deseq_9" , "deseq_10" , "deseq_11" , "deseq_12" , "deseq_13" , "deseq_14" , "deseq_15" , "deseq_16" )
head(cos_1_count)

rv_1_count <- data.frame(sync33[4],sync34[4],sync35[4],sync36[4],sync37[4],sync38[4],sync39[4],sync40[4],sync41[4],sync42[4],sync43[4],sync44[4],sync45[4],sync46[4],sync47[4],sync48[4])
colnames(rv_1_count) <- c( "deseq_33" ,"deseq_34" ,"deseq_35", "deseq_36" , "deseq_37" ,"deseq_38" , "deseq_39" , "deseq_40" , "deseq_41" , "deseq_42" , "deseq_43" , "deseq_44" , "deseq_45" , "deseq_46" , "deseq_47" , "deseq_48" )
head(rv_1_count)


condition = factor (c("0hr" , "3hr" , "6.5hr" , "9hr" , "12hr" , "18.5hr" , "21hr" , "27hr" , "31hr" , "33hr" , "36hr" , "39.5hr" , "42hr" , "45.5hr" , "52hr" , "55hr"))

# cos 1 
sample <- cos_1_count
sample[1:16] <- sapply(sample[1:16], as.integer)
cds = newCountDataSet( sample, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
cos_1_deseq <- counts( cds, normalized=TRUE )



# rv 1 
sample <- rv_1_count
sample[1:16] <- sapply(sample[1:16], as.integer)
cds = newCountDataSet( sample, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
rv_1_deseq <- counts( cds, normalized=TRUE )




##### convert to data frame 
cos_1_deseq <- data.frame(cos_1_deseq)

rv_1_deseq <- data.frame(rv_1_deseq)




sample <- cos_1_deseq
cor_mat <- cor(sample)
png(file="5_figure_correlationmatrix_cos.png",width=10.0,height=6.99,units='in', res=300)
rgb.palette <- colorRampPalette(c("blue", "yellow"), space = "rgb")
levelplot(cor_mat, main="cos1", xlab="", ylab="", col.regions=rgb.palette(120), cuts=100, at=seq(0,1,0.01))
dev.off()


sample <- rv_1_deseq
cor_mat <- cor(sample)
png(file="5_figure_correlationmatrix_rv.png",width=10.0,height=6.99,units='in', res=300)
rgb.palette <- colorRampPalette(c("blue", "yellow"), space = "rgb")
levelplot(cor_mat, main="Rv1", xlab="", ylab="", col.regions=rgb.palette(120), cuts=100, at=seq(0,1,0.01))
dev.off()

sample <- totaldeseq[1:15]
colnames(sample) <- c( "cos2" ,"cos3", "cos4" , "cos5" ,"cos6" , "cos7" , "cos8" , "cos9" , "cos10" , "cos11" , "cos12" , "cos13" , "cos14" , "cos15" , "cos16")
cor_mat <- cor(sample)
png(file="5_figure_correlationmatrix_cos.png",width=10.0,height=6.99,units='in', res=300)
rgb.palette <- colorRampPalette(c("blue", "yellow"), space = "rgb")
levelplot(cor_mat, main="cold sensitive mutant", xlab="", ylab="", col.regions=rgb.palette(120), cuts=100, at=seq(0,1,0.01))
dev.off()


sample <- totaldeseq[31:45]
colnames(sample) <- c( "wt2" ,"wt3", "wt4" , "wt5" ,"wt6" , "wt7" , "wt8" , "wt9" , "wt10" , "wt11" , "wt12" , "wt13" , "wt14" , "wt15" , "wt16")
cor_mat <- cor(sample)
png(file="5_figure_correlationmatrix_rv.png",width=10.0,height=6.99,units='in', res=300)
rgb.palette <- colorRampPalette(c("blue", "yellow"), space = "rgb")
levelplot(cor_mat, main="wild type H37Rv", xlab="", ylab="", col.regions=rgb.palette(120), cuts=100, at=seq(0,1,0.01))
dev.off()




################

cos_2_count <- data.frame(sync17[4],sync18[4],sync19[4],sync20[4],sync21[4],sync22[4],sync23[4],sync24[4],sync25[4],sync26[4],sync27[4],sync28[4],sync29[4],sync30[4],sync31[4],sync32[4])
colnames(cos_2_count) <- c( "deseq_17" ,"deseq_18" ,"deseq_19", "deseq_20" , "deseq_21" ,"deseq_22" , "deseq_23" , "deseq_24" , "deseq_25" , "deseq_26" , "deseq_27" , "deseq_28" , "deseq_29" , "deseq_30" , "deseq_31" , "deseq_32" )
head(cos_2_count)



rv_2_count <- data.frame(sync49[4],sync50[4],sync51[4],sync52[4],sync53[4],sync54[4],sync55[4],sync56[4],sync57[4],sync58[4],sync59[4],sync60[4],sync61[4],sync62[4],sync63[4],sync64[4])
colnames(rv_2_count) <- c( "deseq_49" ,"deseq_50" ,"deseq_51", "deseq_52" , "deseq_53" ,"deseq_54" , "deseq_55" , "deseq_56" , "deseq_57" , "deseq_58" , "deseq_59" , "deseq_60" , "deseq_61" , "deseq_62" , "deseq_63" , "deseq_64" )
head(rv_2_count)


# cos 2 
sample <- cos_2_count
sample[1:16] <- sapply(sample[1:16], as.integer)
cds = newCountDataSet( sample, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
cos_2_deseq <- counts( cds, normalized=TRUE )

# rv 2 
sample <- rv_2_count
sample[1:16] <- sapply(sample[1:16], as.integer)
cds = newCountDataSet( sample, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
rv_2_deseq <- counts( cds, normalized=TRUE )


##### convert to data frame 

cos_2_deseq <- data.frame(cos_2_deseq)
rv_2_deseq <- data.frame(rv_2_deseq)


############################################ cos_1_deseq vs cos_2_deseq scatterplot

p1 <- ggplot(cos_1_deseq, aes(x=cos_1_deseq$deseq_1, y = cos_2_deseq$deseq_17))+geom_point()+xlab("cos1_0")+ylab("cos2_0") + theme(text = element_text(size=20))
p2 <- ggplot(cos_1_deseq, aes(x=cos_1_deseq$deseq_2, y = cos_2_deseq$deseq_18))+geom_point()+xlab("cos1_1")+ylab("cos2_1")+ theme(text = element_text(size=20))
p3 <- ggplot(cos_1_deseq, aes(x=cos_1_deseq$deseq_3, y = cos_2_deseq$deseq_19))+geom_point()+xlab("cos1_2")+ylab("cos2_2")+ theme(text = element_text(size=20))
p4 <- ggplot(cos_1_deseq, aes(x=cos_1_deseq$deseq_4, y = cos_2_deseq$deseq_20))+geom_point()+xlab("cos1_3")+ylab("cos2_3")+ theme(text = element_text(size=20))
p5 <- ggplot(cos_1_deseq, aes(x=cos_1_deseq$deseq_5, y = cos_2_deseq$deseq_21))+geom_point()+xlab("cos1_4")+ylab("cos2_4")+ theme(text = element_text(size=20))
p6 <- ggplot(cos_1_deseq, aes(x=cos_1_deseq$deseq_6, y = cos_2_deseq$deseq_22))+geom_point()+xlab("cos1_5")+ylab("cos2_5")+ theme(text = element_text(size=20))
p7 <- ggplot(cos_1_deseq, aes(x=cos_1_deseq$deseq_7, y = cos_2_deseq$deseq_23))+geom_point()+xlab("cos1_6")+ylab("cos2_6")+ theme(text = element_text(size=20))
p8 <- ggplot(cos_1_deseq, aes(x=cos_1_deseq$deseq_8, y = cos_2_deseq$deseq_24))+geom_point()+xlab("cos1_7")+ylab("cos2_7")+ theme(text = element_text(size=20))
p9 <- ggplot(cos_1_deseq, aes(x=cos_1_deseq$deseq_9, y = cos_2_deseq$deseq_25))+geom_point()+xlab("cos1_8")+ylab("cos2_8")+ theme(text = element_text(size=20))
p10 <- ggplot(cos_1_deseq, aes(x=cos_1_deseq$deseq_10, y = cos_2_deseq$deseq_26))+geom_point()+xlab("cos1_9")+ylab("cos2_9")+ theme(text = element_text(size=20))
p11<- ggplot(cos_1_deseq, aes(x=cos_1_deseq$deseq_11, y = cos_2_deseq$deseq_27))+geom_point()+xlab("cos1_10")+ylab("cos2_10")+ theme(text = element_text(size=20))
p12 <- ggplot(cos_1_deseq, aes(x=cos_1_deseq$deseq_12, y = cos_2_deseq$deseq_28))+geom_point()+xlab("cos1_11")+ylab("cos2_11")+ theme(text = element_text(size=20))
p13 <- ggplot(cos_1_deseq, aes(x=cos_1_deseq$deseq_13, y = cos_2_deseq$deseq_29))+geom_point()+xlab("cos1_12")+ylab("cos2_12")+ theme(text = element_text(size=20))
p14<- ggplot(cos_1_deseq, aes(x=cos_1_deseq$deseq_14, y = cos_2_deseq$deseq_30))+geom_point()+xlab("cos1_13")+ylab("cos2_13")+ theme(text = element_text(size=20))
p15<- ggplot(cos_1_deseq, aes(x=cos_1_deseq$deseq_15, y = cos_2_deseq$deseq_31))+geom_point()+xlab("cos1_14")+ylab("cos2_14")+ theme(text = element_text(size=20))
p16<- ggplot(cos_1_deseq, aes(x=cos_1_deseq$deseq_16, y = cos_2_deseq$deseq_32))+geom_point()+xlab("cos1_15")+ylab("cos2_15")+ theme(text = element_text(size=20))


figure<- grid.arrange(p1,p2,p3,p4, p5,p6,p7,p8,p9,p10,p11,p12, p13,p14,p15,p16, ncol=4, top="")
ggsave("6_figure_cos1_cos2_scatter.png",figure,dpi=300)

p1 <- ggplot(rv_1_deseq, aes(x=rv_1_deseq$deseq_33, y = rv_2_deseq$deseq_49))+geom_point()+xlab("Rv1_0")+ylab("Rv2_0")+ theme(text = element_text(size=20))
p2 <- ggplot(rv_1_deseq, aes(x=rv_1_deseq$deseq_34, y = rv_2_deseq$deseq_50))+geom_point()+xlab("Rv1_1")+ylab("Rv2_1")+ theme(text = element_text(size=20))
p3 <- ggplot(rv_1_deseq, aes(x=rv_1_deseq$deseq_35, y = rv_2_deseq$deseq_51))+geom_point()+xlab("Rv1_2")+ylab("Rv2_2")+ theme(text = element_text(size=20))
p4 <- ggplot(rv_1_deseq, aes(x=rv_1_deseq$deseq_36, y = rv_2_deseq$deseq_52))+geom_point()+xlab("Rv1_3")+ylab("Rv2_3")+ theme(text = element_text(size=20))
p5 <- ggplot(rv_1_deseq, aes(x=rv_1_deseq$deseq_37, y = rv_2_deseq$deseq_53))+geom_point()+xlab("Rv1_4")+ylab("Rv2_4")+ theme(text = element_text(size=20))
p6 <- ggplot(rv_1_deseq, aes(x=rv_1_deseq$deseq_38, y = rv_2_deseq$deseq_54))+geom_point()+xlab("Rv1_5")+ylab("Rv2_5")+ theme(text = element_text(size=20))
p7 <- ggplot(rv_1_deseq, aes(x=rv_1_deseq$deseq_39, y = rv_2_deseq$deseq_55))+geom_point()+xlab("Rv1_6")+ylab("Rv2_6")+ theme(text = element_text(size=20))
p8 <- ggplot(rv_1_deseq, aes(x=rv_1_deseq$deseq_40, y = rv_2_deseq$deseq_56))+geom_point()+xlab("Rv1_7")+ylab("Rv2_7")+ theme(text = element_text(size=20))
p9 <- ggplot(rv_1_deseq, aes(x=rv_1_deseq$deseq_41, y = rv_2_deseq$deseq_57))+geom_point()+xlab("Rv1_8")+ylab("Rv2_8")+ theme(text = element_text(size=20))
p10 <- ggplot(rv_1_deseq, aes(x=rv_1_deseq$deseq_42, y = rv_2_deseq$deseq_58))+geom_point()+xlab("Rv1_9")+ylab("Rv2_9")+ theme(text = element_text(size=20))
p11<- ggplot(rv_1_deseq, aes(x=rv_1_deseq$deseq_43, y = rv_2_deseq$deseq_59))+geom_point()+xlab("Rv1_10")+ylab("Rv2_10")+ theme(text = element_text(size=20))
p12 <- ggplot(rv_1_deseq, aes(x=rv_1_deseq$deseq_44, y = rv_2_deseq$deseq_60))+geom_point()+xlab("Rv1_11")+ylab("Rv2_11")+ theme(text = element_text(size=20))
p13 <- ggplot(rv_1_deseq, aes(x=rv_1_deseq$deseq_45, y = rv_2_deseq$deseq_61))+geom_point()+xlab("Rv1_12")+ylab("Rv2_12")+ theme(text = element_text(size=20))
p14<- ggplot(rv_1_deseq, aes(x=rv_1_deseq$deseq_46, y = rv_2_deseq$deseq_62))+geom_point()+xlab("Rv1_13")+ylab("Rv2_13")+ theme(text = element_text(size=20))
p15<- ggplot(rv_1_deseq, aes(x=rv_1_deseq$deseq_47, y = rv_2_deseq$deseq_63))+geom_point()+xlab("Rv1_14")+ylab("Rv2_14")+ theme(text = element_text(size=20))
p16<- ggplot(rv_1_deseq, aes(x=rv_1_deseq$deseq_48, y = rv_2_deseq$deseq_64))+geom_point()+xlab("Rv1_15")+ylab("Rv2_15")+ theme(text = element_text(size=20))

figure<- grid.arrange(p1,p2,p3,p4, p5,p6,p7,p8,p9,p10,p11,p12, p13,p14,p15,p16, ncol=4, top="")
ggsave("6_figure_rv1_rv2_scatter.png",figure,dpi=300)


#####################################################	global heatmap  ###########################


#dat <- read.csv("7_final_sig_2704_cos_smooth.csv",header=TRUE,row.names=1)
dat <- read.csv("7_rv_smooth_divbymean_adjtp_2704.csv",header=TRUE,row.names=1)

head(dat)

my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
col_breaks = c(seq(0.1,0.6,length=100),  # for red
  seq(0.61,1.0,length=100),           # for yellow
  seq(1.1,1.99,length=100))     

numdata = dat # for cos/rv 2704

colnames(numdata) =c( "3hr" , "6.5hr" , "9hr" , "12hr" , "18.5hr" , "21hr" , "27hr" , "31hr" , "33hr" , "36hr" , "39.5hr" , "42hr" , "45.5hr" , "52hr" , "55hr") 


#png(file="11a_figure_heatmap_siggenes_cos.png",width=7.0,height=7.0,units='in', res=300)
png(file="11b_figure_heatmap_siggenes_rv.png",width=7.0,height=7.0,units='in', res=300)
heatmap.2(as.matrix(numdata), 
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(10,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  Colv="none",
  main ="",
  labRow = FALSE) 


dev.off()

######################################### periodic sin fit #########################
w <- read.csv("sinefit_frequency.csv",header=TRUE)
head(w)
dim(w)
png(file="12_figure_frequencydistribution_2704.png",width=7.0,height=7.0,units='in', res=300)
ggplot(w,aes(w$FREQ))+geom_histogram(bins=100)+xlab("Frequency")+ylab("Distribution")
dev.off()


###############################################################################################
w <- read.csv("fvals.csv",header=TRUE)
head(w)
dim(w)
png(file="13_figure_fvaluedistribution_2165.png",width=7.0,height=7.0,units='in', res=300)
ggplot(w,aes(w$ftest))+geom_histogram(bins=100)+xlab("F-values")+ylab("Distribution")
dev.off()

################## 759 periodic genes clustering ################################



dat <- read.csv("forcluster.csv",header=TRUE,row.names=1)
head(dat)
newdata=dat
disteuc = dist(dat,method="euclidean")
fiteuc = hclust(disteuc,method="ward.D2")
groupseuc = cutree(fiteuc,k=18)
newdata = as.data.frame(newdata)
newdata$clust = groupseuc


colnames(newdata) <- c("gene","3hr" , "6.5hr" , "9hr" , "12hr" , "18.5hr" , "21hr" , "27hr" , "31hr" , "33hr" , "36hr" , "39.5hr" , "42hr" , "45.5hr" , "52hr" , "55hr","clust")



show_graph1 = function(lfcs,cnum) {
  clusttotal=length(rownames(lfcs))
  plot(t(lfcs[1,1:15]),type="l",ylim=c(0,4.0),xlab="",ylab="",main=paste("clust#",cnum,"/ ",clusttotal, "genes"))  
  for (i in 2:length(rownames(lfcs))) {
    lines(t(lfcs[i,1:15])) } 
    med = sapply(lfcs[,1:15],median,2)
    lines(med,col=2) }

png(file="16_figure_cluster_434_periodic.png",width=12.0,height=7.0,units='in', res=300)

par(mfrow=c(3,6))
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


dev.off()
#######################################################################################


############ peak clustering 


dat <- read.csv("onepeak.csv",header = FALSE, row.names=1)

head(dat)

newdata = dat

colnames(newdata) =c("3hr" , "6.5hr" , "9hr" , "12hr" , "18.5hr" , "21hr" , "27hr" , "31hr" , "33hr" , "36hr" , "39.5hr" , "42hr" , "45.5hr" , "52hr" , "55hr", "clust") 



show_graph1 = function(lfcs,cnum) {
  clusttotal=length(rownames(lfcs))
  plot(t(lfcs[1,1:15]),type="l",ylim=c(0,4.0),xlab="",ylab="",main=paste("peak-",cnum,"/ ",clusttotal, "genes"))  
  for (i in 2:length(rownames(lfcs))) {
    lines(t(lfcs[i,1:15])) } 
    med = sapply(lfcs[,1:15],median,2)
    lines(med,col=2) }

png(file="18_figure_clustering_ALLPEAK.png",width=12.0,height=7.0,units='in', res=300)

par(mfrow=c(2,6))
#show_graph1(newdata[newdata$clust==1,],1)
#show_graph1(newdata[newdata$clust==2,],2)
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


dev.off()
##################periodic genes 4 8 10  heatmap 

dat <- read.csv("onepeak.csv",header = FALSE, row.names=1)

head(dat)

newdata = dat





my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
col_breaks = c(seq(0.1,0.6,length=100),  # for red
  seq(0.61,1.0,length=100),           # for yellow
  seq(1.1,1.99,length=100))     


colnames(newdata) =c("3hr" , "6.5hr" , "9hr" , "12hr" , "18.5hr" , "21hr" , "27hr" , "31hr" , "33hr" , "36hr" , "39.5hr" , "42hr" , "45.5hr" , "52hr" , "55hr", "peak") 

#rownames(numdata) <- paste(rownames(dat),dat[,1])

png(file="18_figure_heatmap_PEAK_4.png",width=7.0,height=7.0,units='in', res=300)

heatmap.2(as.matrix(newdata[newdata$peak == 4,1:15]), 
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(10,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  Colv="none",
  main ="") 


dev.off()




png(file="18_figure_heatmap_PEAK_7.png",width=7.0,height=7.0,units='in', res=300)

heatmap.2(as.matrix(newdata[newdata$peak == 7, 1:15]), 
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(10,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  Colv="none",
  main ="") 


dev.off()




png(file="18_figure_heatmap_PEAK_10.png",width=7.0,height=7.0,units='in', res=300)

heatmap.2(as.matrix(newdata[newdata$peak == 10,1:15]), 
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(10,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  Colv="none",
  main ="") 


dev.off()




############ periodic cog category

dat <- read.csv("periodic_cog.csv",header = FALSE)

head(dat)

library(plyr)

count(dat)


dat <- read.csv("significant_cog.csv",header = FALSE)

head(dat)

library(plyr)

count(dat)


#########################################################################################

now final analysis of 2704 significant genes:

clustering ALL SIG GENES:

#dat <- read.csv("7_final_sig_2704_cos_smooth.csv",header=TRUE,row.names=1)
dat <- read.csv("7_rv_smooth_divbymean_adjtp_2704.csv",header=TRUE,row.names=1)


head(dat)
newdata=dat
disteuc = dist(newdata,method="euclidean")
fiteuc = hclust(disteuc,method="ward.D2")
groupseuc = cutree(fiteuc,k=18)
newdata = as.data.frame(newdata)
newdata$clust = groupseuc


colnames(newdata) <- c("3hr" , "6.5hr" , "9hr" , "12hr" , "18.5hr" , "21hr" , "27hr" , "31hr" , "33hr" , "36hr" , "39.5hr" , "42hr" , "45.5hr" , "52hr" , "55hr","clust")



show_graph1 = function(lfcs,cnum) {
  clusttotal=length(rownames(lfcs))
  plot(t(lfcs[1,1:15]),type="l",ylim=c(0,6.0),xlab="",ylab="",main=paste("clust# ",cnum,"/ ",clusttotal, "genes"))  
  for (i in 2:length(rownames(lfcs))) {
    lines(t(lfcs[i,1:15])) } 
    med = sapply(lfcs[,1:15],median,2)
    lines(med,col=2) }

#png(file="20a_figure_cluster_cos_genes_2704_all_sig.png",width=12.0,height=7.0,units='in', res=300)
png(file="20b_figure_cluster_rv_genes_2704_all_sig.png",width=12.0,height=7.0,units='in', res=300)
par(mfrow=c(3,6))
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

dev.off()




dat <- read.csv("6_FINAL_RVProfile_matching_cos_2821.csv",header=TRUE,row.names=1)
head(dat)
newdata=dat
disteuc = dist(newdata,method="euclidean")
fiteuc = hclust(disteuc,method="ward.D2")
groupseuc = cutree(fiteuc,k=18)
newdata = as.data.frame(newdata)
newdata$clust = groupseuc


colnames(newdata) <- c("3hr" , "6.5hr" , "9hr" , "12hr" , "18.5hr" , "21hr" , "27hr" , "31hr" , "33hr" , "36hr" , "39.5hr" , "42hr" , "45.5hr" , "52hr" , "55hr","clust")



show_graph1 = function(lfcs,cnum) {
  clusttotal=length(rownames(lfcs))
  plot(t(lfcs[1,1:15]),type="l",ylim=c(0,6.0),xlab="",ylab="",main=paste("clust# ",cnum,"/ ",clusttotal, "genes"))  
  for (i in 2:length(rownames(lfcs))) {
    lines(t(lfcs[i,1:15])) } 
    med = sapply(lfcs[,1:15],median,2)
    lines(med,col=2) }

png(file="20b_figure_cluster_genes_2821_all_sig_rv.png",width=12.0,height=7.0,units='in', res=300)

par(mfrow=c(3,6))
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

dev.off()









#########################################################################################

############# 

#####################################

CLUSTER TFS:

dat <- read.csv("test.csv",header=TRUE,row.names=1)
head(dat)
newdata=dat[1:16]
disteuc = dist(newdata,method="euclidean")
fiteuc = hclust(disteuc,method="ward.D2")
groupseuc = cutree(fiteuc,k=6)
newdata = as.data.frame(newdata)
newdata$clust = groupseuc


colnames(newdata) <- c("3hr" , "6.5hr" , "9hr" , "12hr" , "18.5hr" , "21hr" , "27hr" , "31hr" , "33hr" , "36hr" , "39.5hr" , "42hr" , "45.5hr" , "52hr" , "55hr","clust")



show_graph1 = function(lfcs,cnum) {
  clusttotal=length(rownames(lfcs))
  plot(t(lfcs[1,1:15]),type="l",ylim=c(0,6.0),xlab="",ylab="",main=paste("clust# ",cnum,"/ ",clusttotal, "genes"))  
  for (i in 2:length(rownames(lfcs))) {
    lines(t(lfcs[i,1:15])) } 
    med = sapply(lfcs[,1:15],median,2)
    lines(med,col=2) }

#png(file="20_figure_cluster_genes_617.png",width=12.0,height=7.0,units='in', res=300)

par(mfrow=c(2,6))
show_graph1(newdata[newdata$clust==0,],0)
#show_graph1(newdata[newdata$clust==2,],2)
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
#show_graph1(newdata[newdata$clust==14,],14)
#show_graph1(newdata[newdata$clust==15,],15)
#show_graph1(newdata[newdata$clust==16,],16)
#show_graph1(newdata[newdata$clust==17,],17)
#show_graph1(newdata[newdata$clust==18,],18)

#dev.off()

dat <- read.csv("1_TF_profiles_forCLUSTERING.csv",header=FALSE,row.names=1)
head(dat)

numdata = dat[,2:16] 


my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
col_breaks = c(seq(0.1,0.6,length=100),  # for red
  seq(0.61,1.0,length=100),           # for yellow
  seq(1.1,1.99,length=100))     


colnames(numdata) =c( "3hr" , "6.5hr" , "9hr" , "12hr" , "18.5hr" , "21hr" , "27hr" , "31hr" , "33hr" , "36hr" , "39.5hr" , "42hr" , "45.5hr" , "52hr" , "55hr") 

rownames(numdata)= paste(rownames(dat),dat[,1])

png(file="222_figure_heatmap_tfs.png",width=10.0,height=10.0,units='in', res=300)

heatmap.2(as.matrix(numdata), 
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(10,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  Colv="none",
  main ="") 
dev.off()


# Install factoextra
install.packages("factoextra")
# Install cluster package
install.packages("cluster")

library("cluster")
library("factoextra")

head(newdata)

#find the optimal number of cluster
fviz_nbclust(newdata, kmeans, method = "gap_stat")



############### Compute PAM   ################33
library("cluster")
pam.res <- pam(newdata, 5)

# Visualize
fviz_cluster(pam.res)


############### k means ##################
km.res <- kmeans(newdata[1:15], 10, nstart = 25)
# Visualize
library("factoextra")
fviz_cluster(km.res, data = newdata, frame.type = "convex")+
  theme_minimal()

