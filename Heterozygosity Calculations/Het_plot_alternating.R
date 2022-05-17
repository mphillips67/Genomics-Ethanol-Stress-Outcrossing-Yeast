
library(GenWin)
library(dplyr)
par(mfrow=c(6,1),mar=c(2,5,0.5,2),mgp=c(3,1,0), oma=c(3,4,0,0))
#treatment + tps
control <-c(paste("0",1:9,sep=""),10,12:21)
moderate <- c(25:30,32:33,35:38,41:48)
high <- c(49:56,59:61,64:72)

#snp table
data <- read.table("SNPtable_50X_scaled.txt", header = TRUE)
maf <- data[,c(seq(5,ncol(data),by=2))]
hets <- maf
for (i in 1:ncol(maf)){
  hets[,i] <- 2 * hets[,i] * (1-hets[,i])
}
hets <- cbind(data[,1:2],hets)


#control het week 1
control_week15 <- paste("alt_ETH_rep", control,"_wk01",sep="")
het_cont_week15 <- cbind(hets[,1:2],rowMeans(hets[,control_week15]))
names(het_cont_week15) <- c("chr","pos","het")
Control_wk15_final <- data.frame()
chr_index <- unique(het_cont_week15$chr)
for (i in 1:16){
  
  temp <- subset(het_cont_week15, chr == chr_index[i])
  
  temp2 <- splineAnalyze(temp$het,temp$pos,smoothness=3000,
                         plotRaw=FALSE,plotWindows=FALSE,method=3 )
  
  het <- temp2[["windowData"]][["MeanY"]]
  pos <- temp2[["windowData"]][["WindowStop"]]
  chr <- rep( chr_index[i], length(het))
  temp3 <- data.frame(chr,pos,het )
  Control_wk15_final <- rbind(Control_wk15_final, temp3)
}


#Make a genome-wide scale so the entire genome can go on a single plot
#Add the max pos of C01 to C02, C01+C02 to C03, etc
#Asssuming you are starting with a data frame called xx with one colname = pos, and another colname = var1 (some variable you want to plot)
Gaxis <- numeric(length=0)
chrs<-c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')

Control_wk15_final <- arrange(transform(Control_wk15_final,
                                        chr=factor(chr,levels=chrs)),chr)

for (k in 1:length(chrs)){
  data.samp<-subset(Control_wk15_final,chr==chrs[k])
  if (chrs[k]=='C01'){Gaxis.samp<-data.samp$pos}
  if (chrs[k]=='C02'){Gaxis.samp<-data.samp$pos+230218}
  if (chrs[k]=='C03'){Gaxis.samp<-data.samp$pos+230218+813184}
  if (chrs[k]=='C04'){Gaxis.samp<-data.samp$pos+230218+813184+316620}
  if (chrs[k]=='C05'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933}
  if (chrs[k]=='C06'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874}
  if (chrs[k]=='C07'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161}
  if (chrs[k]=='C08'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940}
  if (chrs[k]=='C09'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643}
  if (chrs[k]=='C10'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888}
  if (chrs[k]=='C11'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751}
  if (chrs[k]=='C12'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816} 
  if (chrs[k]=='C13'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177}
  if (chrs[k]=='C14'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431}
  if (chrs[k]=='C15'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333}
  if (chrs[k]=='C16'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333+1091291}
  Gaxis<-c(Gaxis,Gaxis.samp) 
}
MB=Gaxis/1e6

plot(MB,Control_wk15_final$het,
     xlab="",
     ylab="",
     main="",
     ylim=c(0,0.5),
     cex.main=3, #size of text for main title
     cex.lab=2, #size of text for axis labels
     type="n",	#type = none: don't make plot yet, need to make rectangles first
     axes=FALSE  #we will do the axes by hand next
)

#Here draw are gray rectangles that delineate chromosomes
rect(0.23,-10,1.04,80,col="grey80",lty=0)  
rect(1.36,-10,2.89,80,col="grey80",lty=0)  
rect(3.46,-10,3.73,80,col="grey80",lty=0)
rect(4.82,-10,5.39,80,col="grey80",lty=0)  
rect(5.83,-10,6.57,80,col="grey80",lty=0) 
rect(7.24,-10,8.32,80,col="grey80",lty=0)
rect(9.24,-10,10.03,80,col="grey80",lty=0)
rect(11.12,-10,12.07,80,col="grey80",lty=0) 

#Here draw points or lines.  Often it is useful to play around with point type, font size, color...
points(MB,Control_wk15_final$het, col="darkgreen", pch = 20, cex = 0.60, type = "l")
box() 
#Now draw axes back in (you have more flexibility this way)
axis(2,cex.axis=1.5, las = 1)

#legend("right",inset=-.07, "C",bty="n", cex=2)
legend("left",inset=-.1, "(A)",bty="n", cex=2)


#control cycle 15
control_week15 <- paste("alt_ETH_rep", control,"_wk15",sep="")
het_cont_week15 <- cbind(hets[,1:2],rowMeans(hets[,control_week15]))
names(het_cont_week15) <- c("chr","pos","het")
Control_wk15_final <- data.frame()
chr_index <- unique(het_cont_week15$chr)
for (i in 1:16){
  
  temp <- subset(het_cont_week15, chr == chr_index[i])
  
  temp2 <- splineAnalyze(temp$het,temp$pos,smoothness=3000,
                         plotRaw=FALSE,plotWindows=FALSE,method=3 )
  
  het <- temp2[["windowData"]][["MeanY"]]
  pos <- temp2[["windowData"]][["WindowStop"]]
  chr <- rep( chr_index[i], length(het))
  temp3 <- data.frame(chr,pos,het )
  Control_wk15_final <- rbind(Control_wk15_final, temp3)
}


#Make a genome-wide scale so the entire genome can go on a single plot
#Add the max pos of C01 to C02, C01+C02 to C03, etc
#Asssuming you are starting with a data frame called xx with one colname = pos, and another colname = var1 (some variable you want to plot)
Gaxis <- numeric(length=0)
chrs<-c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')

Control_wk15_final <- arrange(transform(Control_wk15_final,
                                        chr=factor(chr,levels=chrs)),chr)

for (k in 1:length(chrs)){
  data.samp<-subset(Control_wk15_final,chr==chrs[k])
  if (chrs[k]=='C01'){Gaxis.samp<-data.samp$pos}
  if (chrs[k]=='C02'){Gaxis.samp<-data.samp$pos+230218}
  if (chrs[k]=='C03'){Gaxis.samp<-data.samp$pos+230218+813184}
  if (chrs[k]=='C04'){Gaxis.samp<-data.samp$pos+230218+813184+316620}
  if (chrs[k]=='C05'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933}
  if (chrs[k]=='C06'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874}
  if (chrs[k]=='C07'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161}
  if (chrs[k]=='C08'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940}
  if (chrs[k]=='C09'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643}
  if (chrs[k]=='C10'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888}
  if (chrs[k]=='C11'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751}
  if (chrs[k]=='C12'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816} 
  if (chrs[k]=='C13'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177}
  if (chrs[k]=='C14'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431}
  if (chrs[k]=='C15'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333}
  if (chrs[k]=='C16'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333+1091291}
  Gaxis<-c(Gaxis,Gaxis.samp) 
}
MB=Gaxis/1e6

plot(MB,Control_wk15_final$het,
     xlab="",
     ylab="",
     main="",
     ylim=c(0,0.5),
     cex.main=3, #size of text for main title
     cex.lab=2, #size of text for axis labels
     type="n",	#type = none: don't make plot yet, need to make rectangles first
     axes=FALSE  #we will do the axes by hand next
)

#Here draw are gray rectangles that delineate chromosomes
rect(0.23,-10,1.04,80,col="grey80",lty=0)  
rect(1.36,-10,2.89,80,col="grey80",lty=0)  
rect(3.46,-10,3.73,80,col="grey80",lty=0)
rect(4.82,-10,5.39,80,col="grey80",lty=0)  
rect(5.83,-10,6.57,80,col="grey80",lty=0) 
rect(7.24,-10,8.32,80,col="grey80",lty=0)
rect(9.24,-10,10.03,80,col="grey80",lty=0)
rect(11.12,-10,12.07,80,col="grey80",lty=0) 

#Here draw points or lines.  Often it is useful to play around with point type, font size, color...
points(MB,Control_wk15_final$het, col="darkgreen", pch = 20, cex = 0.60, type = "l")
box() 
#Now draw axes back in (you have more flexibility this way)
axis(2,cex.axis=1.5, las = 1)

#legend("right",inset=-.07, "C",bty="n", cex=2)
legend("left",inset=-.1, "(B)",bty="n", cex=2)


#moderarate c1
moderate_week15 <- paste("alt_ETH_rep", moderate,"_wk01",sep="")
het_mod_week15 <- cbind(hets[,1:2],rowMeans(hets[,moderate_week15]))
names(het_mod_week15) <- c("chr","pos","het")
moderate_wk15_final <- data.frame()
chr_index <- unique(het_mod_week15$chr)
for (i in 1:16){
  
  temp <- subset(het_mod_week15, chr == chr_index[i])
  
  temp2 <- splineAnalyze(temp$het,temp$pos,smoothness=3000,
                         plotRaw=FALSE,plotWindows=FALSE,method=3 )
  
  het <- temp2[["windowData"]][["MeanY"]]
  pos <- temp2[["windowData"]][["WindowStop"]]
  chr <- rep( chr_index[i], length(het))
  temp3 <- data.frame(chr,pos,het )
  moderate_wk15_final <- rbind(moderate_wk15_final, temp3)
}


#Make a genome-wide scale so the entire genome can go on a single plot
#Add the max pos of C01 to C02, C01+C02 to C03, etc
#Asssuming you are starting with a data frame called xx with one colname = pos, and another colname = var1 (some variable you want to plot)
Gaxis <- numeric(length=0)
chrs<-c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')

moderate_wk15_final <- arrange(transform(moderate_wk15_final,
                                         chr=factor(chr,levels=chrs)),chr)

for (k in 1:length(chrs)){
  data.samp<-subset(moderate_wk15_final,chr==chrs[k])
  if (chrs[k]=='C01'){Gaxis.samp<-data.samp$pos}
  if (chrs[k]=='C02'){Gaxis.samp<-data.samp$pos+230218}
  if (chrs[k]=='C03'){Gaxis.samp<-data.samp$pos+230218+813184}
  if (chrs[k]=='C04'){Gaxis.samp<-data.samp$pos+230218+813184+316620}
  if (chrs[k]=='C05'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933}
  if (chrs[k]=='C06'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874}
  if (chrs[k]=='C07'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161}
  if (chrs[k]=='C08'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940}
  if (chrs[k]=='C09'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643}
  if (chrs[k]=='C10'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888}
  if (chrs[k]=='C11'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751}
  if (chrs[k]=='C12'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816} 
  if (chrs[k]=='C13'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177}
  if (chrs[k]=='C14'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431}
  if (chrs[k]=='C15'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333}
  if (chrs[k]=='C16'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333+1091291}
  Gaxis<-c(Gaxis,Gaxis.samp) 
}
MB=Gaxis/1e6

plot(MB,moderate_wk15_final$het,
     xlab="",
     ylab="",
     main="",
     ylim=c(0,0.5),
     cex.main=3, #size of text for main title
     cex.lab=2, #size of text for axis labels
     type="n",	#type = none: don't make plot yet, need to make rectangles first
     axes=FALSE  #we will do the axes by hand next
)

#Here draw are gray rectangles that delineate chromosomes
rect(0.23,-10,1.04,80,col="grey80",lty=0)  
rect(1.36,-10,2.89,80,col="grey80",lty=0)  
rect(3.46,-10,3.73,80,col="grey80",lty=0)
rect(4.82,-10,5.39,80,col="grey80",lty=0)  
rect(5.83,-10,6.57,80,col="grey80",lty=0) 
rect(7.24,-10,8.32,80,col="grey80",lty=0)
rect(9.24,-10,10.03,80,col="grey80",lty=0)
rect(11.12,-10,12.07,80,col="grey80",lty=0) 

#Here draw points or lines.  Often it is useful to play around with point type, font size, color...
points(MB,moderate_wk15_final$het, col="darkorange", pch = 20, cex = 0.60, type = "l")
box() 
#Now draw axes back in (you have more flexibility this way)
axis(2,cex.axis=1.5, las = 1)

#legend("right",inset=-.07, "C",bty="n", cex=2)
legend("left",inset=-.1, "(C)",bty="n", cex=2)


#moderate het week 15
moderate_week15 <- paste("alt_ETH_rep", moderate,"_wk15",sep="")
het_mod_week15 <- cbind(hets[,1:2],rowMeans(hets[,moderate_week15]))
names(het_mod_week15) <- c("chr","pos","het")
moderate_wk15_final <- data.frame()
chr_index <- unique(het_mod_week15$chr)
for (i in 1:16){
  
  temp <- subset(het_mod_week15, chr == chr_index[i])
  
  temp2 <- splineAnalyze(temp$het,temp$pos,smoothness=3000,
                         plotRaw=FALSE,plotWindows=FALSE,method=3 )
  
  het <- temp2[["windowData"]][["MeanY"]]
  pos <- temp2[["windowData"]][["WindowStop"]]
  chr <- rep( chr_index[i], length(het))
  temp3 <- data.frame(chr,pos,het )
  moderate_wk15_final <- rbind(moderate_wk15_final, temp3)
}


#Make a genome-wide scale so the entire genome can go on a single plot
#Add the max pos of C01 to C02, C01+C02 to C03, etc
#Asssuming you are starting with a data frame called xx with one colname = pos, and another colname = var1 (some variable you want to plot)
Gaxis <- numeric(length=0)
chrs<-c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')

moderate_wk15_final <- arrange(transform(moderate_wk15_final,
                                         chr=factor(chr,levels=chrs)),chr)

for (k in 1:length(chrs)){
  data.samp<-subset(moderate_wk15_final,chr==chrs[k])
  if (chrs[k]=='C01'){Gaxis.samp<-data.samp$pos}
  if (chrs[k]=='C02'){Gaxis.samp<-data.samp$pos+230218}
  if (chrs[k]=='C03'){Gaxis.samp<-data.samp$pos+230218+813184}
  if (chrs[k]=='C04'){Gaxis.samp<-data.samp$pos+230218+813184+316620}
  if (chrs[k]=='C05'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933}
  if (chrs[k]=='C06'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874}
  if (chrs[k]=='C07'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161}
  if (chrs[k]=='C08'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940}
  if (chrs[k]=='C09'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643}
  if (chrs[k]=='C10'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888}
  if (chrs[k]=='C11'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751}
  if (chrs[k]=='C12'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816} 
  if (chrs[k]=='C13'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177}
  if (chrs[k]=='C14'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431}
  if (chrs[k]=='C15'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333}
  if (chrs[k]=='C16'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333+1091291}
  Gaxis<-c(Gaxis,Gaxis.samp) 
}
MB=Gaxis/1e6

plot(MB,moderate_wk15_final$het,
     xlab="",
     ylab="",
     main="",
     ylim=c(0,0.5),
     cex.main=3, #size of text for main title
     cex.lab=2, #size of text for axis labels
     type="n",	#type = none: don't make plot yet, need to make rectangles first
     axes=FALSE  #we will do the axes by hand next
)

#Here draw are gray rectangles that delineate chromosomes
rect(0.23,-10,1.04,80,col="grey80",lty=0)  
rect(1.36,-10,2.89,80,col="grey80",lty=0)  
rect(3.46,-10,3.73,80,col="grey80",lty=0)
rect(4.82,-10,5.39,80,col="grey80",lty=0)  
rect(5.83,-10,6.57,80,col="grey80",lty=0) 
rect(7.24,-10,8.32,80,col="grey80",lty=0)
rect(9.24,-10,10.03,80,col="grey80",lty=0)
rect(11.12,-10,12.07,80,col="grey80",lty=0) 

#Here draw points or lines.  Often it is useful to play around with point type, font size, color...
points(MB,moderate_wk15_final$het, col="darkorange", pch = 20, cex = 0.60, type = "l")
box() 
#Now draw axes back in (you have more flexibility this way)
axis(2,cex.axis=1.5, las = 1)

#legend("right",inset=-.07, "C",bty="n", cex=2)
legend("left",inset=-.1, "(D)",bty="n", cex=2)




#high het week 15
high_week15 <- paste("alt_ETH_rep", high,"_wk01",sep="")
het_high_week15 <- cbind(hets[,1:2],rowMeans(hets[,high_week15]))
names(het_high_week15) <- c("chr","pos","het")
high_wk15_final <- data.frame()
chr_index <- unique(het_high_week15$chr)
for (i in 1:16){
  
  temp <- subset(het_high_week15, chr == chr_index[i])
  
  temp2 <- splineAnalyze(temp$het,temp$pos,smoothness=3000,
                         plotRaw=FALSE,plotWindows=FALSE,method=3 )
  
  het <- temp2[["windowData"]][["MeanY"]]
  pos <- temp2[["windowData"]][["WindowStop"]]
  chr <- rep( chr_index[i], length(het))
  temp3 <- data.frame(chr,pos,het )
  high_wk15_final <- rbind(high_wk15_final, temp3)
}


#Make a genome-wide scale so the entire genome can go on a single plot
#Add the max pos of C01 to C02, C01+C02 to C03, etc
#Asssuming you are starting with a data frame called xx with one colname = pos, and another colname = var1 (some variable you want to plot)
Gaxis <- numeric(length=0)
chrs<-c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')

high_wk15_final <- arrange(transform(high_wk15_final,
                                     chr=factor(chr,levels=chrs)),chr)

for (k in 1:length(chrs)){
  data.samp<-subset(high_wk15_final,chr==chrs[k])
  if (chrs[k]=='C01'){Gaxis.samp<-data.samp$pos}
  if (chrs[k]=='C02'){Gaxis.samp<-data.samp$pos+230218}
  if (chrs[k]=='C03'){Gaxis.samp<-data.samp$pos+230218+813184}
  if (chrs[k]=='C04'){Gaxis.samp<-data.samp$pos+230218+813184+316620}
  if (chrs[k]=='C05'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933}
  if (chrs[k]=='C06'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874}
  if (chrs[k]=='C07'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161}
  if (chrs[k]=='C08'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940}
  if (chrs[k]=='C09'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643}
  if (chrs[k]=='C10'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888}
  if (chrs[k]=='C11'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751}
  if (chrs[k]=='C12'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816} 
  if (chrs[k]=='C13'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177}
  if (chrs[k]=='C14'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431}
  if (chrs[k]=='C15'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333}
  if (chrs[k]=='C16'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333+1091291}
  Gaxis<-c(Gaxis,Gaxis.samp) 
}
MB=Gaxis/1e6

plot(MB,high_wk15_final$het,
     xlab="",
     ylab="",
     main="",
     ylim=c(0,0.5),
     cex.main=3, #size of text for main title
     cex.lab=2, #size of text for axis labels
     type="n",	#type = none: don't make plot yet, need to make rectangles first
     axes=FALSE  #we will do the axes by hand next
)

#Here draw are gray rectangles that delineate chromosomes
rect(0.23,-10,1.04,80,col="grey80",lty=0)  
rect(1.36,-10,2.89,80,col="grey80",lty=0)  
rect(3.46,-10,3.73,80,col="grey80",lty=0)
rect(4.82,-10,5.39,80,col="grey80",lty=0)  
rect(5.83,-10,6.57,80,col="grey80",lty=0) 
rect(7.24,-10,8.32,80,col="grey80",lty=0)
rect(9.24,-10,10.03,80,col="grey80",lty=0)
rect(11.12,-10,12.07,80,col="grey80",lty=0) 

#Here draw points or lines.  Often it is useful to play around with point type, font size, color...
points(MB,high_wk15_final$het, col="darkred", pch = 20, cex = 0.60, type = "l")
box() 
#Now draw axes back in (you have more flexibility this way)
axis(2,cex.axis=1.5, las = 1)

#legend("right",inset=-.07, "C",bty="n", cex=2)
legend("left",inset=-.1, "(E)",bty="n", cex=2)



#high het week 15
high_week15 <- paste("alt_ETH_rep", high,"_wk15",sep="")
het_high_week15 <- cbind(hets[,1:2],rowMeans(hets[,high_week15]))
names(het_high_week15) <- c("chr","pos","het")
high_wk15_final <- data.frame()
chr_index <- unique(het_high_week15$chr)
for (i in 1:16){
  
  temp <- subset(het_high_week15, chr == chr_index[i])
  
  temp2 <- splineAnalyze(temp$het,temp$pos,smoothness=3000,
                         plotRaw=FALSE,plotWindows=FALSE,method=3 )
  
  het <- temp2[["windowData"]][["MeanY"]]
  pos <- temp2[["windowData"]][["WindowStop"]]
  chr <- rep( chr_index[i], length(het))
  temp3 <- data.frame(chr,pos,het )
  high_wk15_final <- rbind(high_wk15_final, temp3)
}


#Make a genome-wide scale so the entire genome can go on a single plot
#Add the max pos of C01 to C02, C01+C02 to C03, etc
#Asssuming you are starting with a data frame called xx with one colname = pos, and another colname = var1 (some variable you want to plot)
Gaxis <- numeric(length=0)
chrs<-c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')

high_wk15_final <- arrange(transform(high_wk15_final,
                                     chr=factor(chr,levels=chrs)),chr)

for (k in 1:length(chrs)){
  data.samp<-subset(high_wk15_final,chr==chrs[k])
  if (chrs[k]=='C01'){Gaxis.samp<-data.samp$pos}
  if (chrs[k]=='C02'){Gaxis.samp<-data.samp$pos+230218}
  if (chrs[k]=='C03'){Gaxis.samp<-data.samp$pos+230218+813184}
  if (chrs[k]=='C04'){Gaxis.samp<-data.samp$pos+230218+813184+316620}
  if (chrs[k]=='C05'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933}
  if (chrs[k]=='C06'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874}
  if (chrs[k]=='C07'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161}
  if (chrs[k]=='C08'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940}
  if (chrs[k]=='C09'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643}
  if (chrs[k]=='C10'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888}
  if (chrs[k]=='C11'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751}
  if (chrs[k]=='C12'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816} 
  if (chrs[k]=='C13'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177}
  if (chrs[k]=='C14'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431}
  if (chrs[k]=='C15'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333}
  if (chrs[k]=='C16'){Gaxis.samp<-data.samp$pos+230218+813184+316620+1531933+576874+270161+1090940+562643+439888+745751+666816+1078177+924431+784333+1091291}
  Gaxis<-c(Gaxis,Gaxis.samp) 
}
MB=Gaxis/1e6

plot(MB,high_wk15_final$het,
     xlab="",
     ylab="",
     main="",
     ylim=c(0,0.5),
     cex.main=3, #size of text for main title
     cex.lab=2, #size of text for axis labels
     type="n",	#type = none: don't make plot yet, need to make rectangles first
     axes=FALSE  #we will do the axes by hand next
)

#Here draw are gray rectangles that delineate chromosomes
rect(0.23,-10,1.04,80,col="grey80",lty=0)  
rect(1.36,-10,2.89,80,col="grey80",lty=0)  
rect(3.46,-10,3.73,80,col="grey80",lty=0)
rect(4.82,-10,5.39,80,col="grey80",lty=0)  
rect(5.83,-10,6.57,80,col="grey80",lty=0) 
rect(7.24,-10,8.32,80,col="grey80",lty=0)
rect(9.24,-10,10.03,80,col="grey80",lty=0)
rect(11.12,-10,12.07,80,col="grey80",lty=0) 

#Here draw points or lines.  Often it is useful to play around with point type, font size, color...
points(MB,high_wk15_final$het, col="darkred", pch = 20, cex = 0.60, type = "l")
box() 
#Now draw axes back in (you have more flexibility this way)
axis(2,cex.axis=1.5, las = 1)

#legend("right",inset=-.07, "C15",bty="n", cex=2)
legend("left",inset=-.1, "(F)",bty="n", cex=2)

#add y and x axis label
mtext(text="Heterozygosity",side=2,line=0,outer=TRUE,cex=1.5)
mtext(text="Position (mb)",side=1,line=1.5,outer=TRUE,cex=1.5)

#add chromosome labels at midpts
mtext("C2",line = .5,side=1, at =0.635, cex=1.2)
mtext("C4",line = .5,side=1, at =2.125, cex=1.2)
mtext("C6",line = .5,side=1, at =3.595, cex=1.2)
mtext("C8",line = .5,side=1, at =5.105, cex=1.2)
mtext("C10",line = .5,side=1, at =6.2, cex=1.2)
mtext("C12",line = .5,side=1, at =7.78, cex=1.2)
mtext("C14",line = .5,side=1, at =9.635, cex=1.2)
mtext("C16",line = .5,side=1, at =11.595, cex=1.2)



