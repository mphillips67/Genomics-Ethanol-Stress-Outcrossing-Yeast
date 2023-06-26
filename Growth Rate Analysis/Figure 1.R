library(ggplot2)
library(reshape2)
library(gridExtra)
library("matrixStats")
data <- read.table("Plain_ypd_2021.txt", header = TRUE)
data[,2:ncol(data)] <- log10(data[,2:ncol(data)])
data <- data[2:nrow(data),]
ancestor <- data[,3:11]
control <- data[,12:21]
e6 <- data[,22:31]
e10 <- data[,32:41]
sd <- c(rowSds(as.matrix(ancestor)),rowSds(as.matrix(control)),rowSds(as.matrix(e6)),rowSds(as.matrix(e10)))
means <- data.frame(data$time/60,rowMeans(ancestor),rowMeans(control),rowMeans(e6),rowMeans(e10))
names(means) <- c("time", "Ancestor", "C Replicates", "M Replicates","H Replicates")
data2 <- melt(means, id.vars="time")
data2 <- cbind(data2, sd)
names(data2) <- c("Time","Treatment","OD","sd")
data2$Treatment <- factor(data2$Treatment, levels = c("Ancestor", "C Replicates", "M Replicates","H Replicates"))
a <- ggplot(data2, aes(x=Time, y=OD,col=Treatment)) +  ylim(-1,.25) +
  geom_point() +  geom_errorbar( aes(ymin = OD-sd, ymax = OD+sd),width = 0.2) + ylab('log(OD)') +  xlab('Time (hr)')  + scale_color_manual(breaks=c("Ancestor", "C Replicates", "M Replicates","H Replicates"), values =c("black","#4daf4a","#ff7f00","#e41a1c"))  + ggtitle("A. Growth Curves in Plain YPD") + theme( axis.text=element_text(size=10),
                                                                                                                                                                                                                                                                                                        axis.title=element_text(size=10,face="bold")) + theme(legend.title = element_text(size=10, face="bold")) + theme(legend.position="bottom") 
#doubling time 
library(growthcurver)
library(growthcurver)
data <- read.table("Plain_ypd_2021.txt", header = TRUE)
data <- data[2:nrow(data),]
data$time <- data$time/60
fit_all <- SummarizeGrowthByPlate(data, bg_correct = "blank")
anc <- c(fit_all[1:9,6], NA)
doubling_time <-  data.frame(anc,fit_all[10:19,6],fit_all[20:29,6],fit_all[30:39,6])
Treatment <- c(rep("Ancestor", nrow(doubling_time)),rep("C Replicates", nrow(doubling_time)),rep("M Replicates", nrow(doubling_time)),rep("H Replicates", nrow(doubling_time)))
DD <-   c(anc,fit_all[10:19,6],fit_all[20:29,6],fit_all[30:39,6])
data3 <- data.frame(Treatment, DD)
data3$Treatment <- factor(data3$Treatment, levels = c("Ancestor", "C Replicates", "M Replicates","H Replicates"))
b <- ggplot(data3, aes(x=Treatment, y=DD, fill=Treatment)) +
  ylab('Doubling Time (hr)') +  ggtitle("B. Doubling Time in Plain YPD") +
  geom_boxplot()  + scale_fill_manual(breaks=c("Ancestor", "C Replicates", "M Replicates","H Replicates"), values =c("black", "#4daf4a","#ff7f00","#e41a1c")) + ylim(1,1.5) + theme( axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold")) + theme(legend.title = element_text(size=15, face="bold"))  + theme(legend.position="none")


#moderate
library(ggplot2)
library(reshape2)
library("matrixStats")
data <- read.table("Moderate_ethanol_2021.txt", header = TRUE)
data[,2:ncol(data)] <- log10(data[,2:ncol(data)])
data <- data[2:nrow(data),]

ancestor <- data[,3:11]
control <- data[,12:21]
e6 <- data[,22:31]
e10 <- data[,32:41]
sd <- c(rowSds(as.matrix(ancestor)),rowSds(as.matrix(control)),rowSds(as.matrix(e6)),rowSds(as.matrix(e10)))

means <- data.frame(data$time/60,rowMeans(ancestor),rowMeans(control),rowMeans(e6),rowMeans(e10))
names(means) <- c("time", "Ancestor", "C Replicates", "M Replicates","H Replicates")
data2 <- melt(means, id.vars="time")
data2 <- cbind(data2, sd)

names(data2) <- c("Time","Treatment","OD","sd")
data2$Treatment <- factor(data2$Treatment, levels = c("Ancestor", "C Replicates", "M Replicates","H Replicates"))

c <- ggplot(data2, aes(x=Time, y=OD,col=Treatment)) + ylim(-1,.3) +
  geom_point() +  geom_errorbar( aes(ymin = OD-sd, ymax = OD+sd),width = 0.25) + ylab('log(OD)') +  xlab('Time (hr)')  + scale_color_manual(breaks=c("Ancestor", "C Replicates", "M Replicates","H Replicates"), values =c("black","#4daf4a","#ff7f00","#e41a1c"))  + ggtitle("C. Growth Curves in 6% Ethanol YPD") + theme( axis.text=element_text(size=10),
                                                                                                                                                                                                                                                                                                                             axis.title=element_text(size=10,face="bold")) + theme(legend.title = element_text(size=10, face="bold")) + theme(legend.position="bottom") 

#doubling time 
library(growthcurver)
library(growthcurver)
data <- read.table("Moderate_ethanol_2021.txt", header = TRUE)
data <- data[2:nrow(data),]
data$time <- data$time/60
fit_all <- SummarizeGrowthByPlate(data, bg_correct = "blank")
anc <- c(fit_all[1:9,6], NA)
doubling_time <-  data.frame(anc,fit_all[10:19,6],fit_all[20:29,6],fit_all[30:39,6])
Treatment <- c(rep("Ancestor", nrow(doubling_time)),rep("C Replicates", nrow(doubling_time)),rep("M Replicates", nrow(doubling_time)),rep("H Replicates", nrow(doubling_time)))
DD <-   c(anc,fit_all[10:19,6],fit_all[20:29,6],fit_all[30:39,6])
data3 <- data.frame(Treatment, DD)


data3$Treatment <- factor(data3$Treatment, levels = c("Ancestor", "C Replicates", "M Replicates","H Replicates"))
d <- ggplot(data3, aes(x=Treatment, y=DD, fill=Treatment)) +
  ylab('Doubling Time (hr)') +  ggtitle("D. Doubling Time in 6% Ethanol YPD") +
  geom_boxplot()  + scale_fill_manual(breaks=c("Ancestor", "C Replicates", "M Replicates","H Replicates"), values =c("black", "#4daf4a","#ff7f00","#e41a1c")) + ylim(1,2.5) + theme( axis.text=element_text(size=10),
                                                                                                                                                                                   axis.title=element_text(size=10,face="bold")) + theme(legend.title = element_text(size=15, face="bold")) + theme(legend.position="none")

#high
library(ggplot2)
library(reshape2)
library("matrixStats")
data <- read.table("High_ethanol_2021.txt", header = TRUE)
data[,2:ncol(data)] <- log10(data[,2:ncol(data)])
data <- data[2:nrow(data),]

ancestor <- data[,3:11]
control <- data[,12:21]
e6 <- data[,22:31]
e10 <- data[,32:41]
sd <- c(rowSds(as.matrix(ancestor)),rowSds(as.matrix(control)),rowSds(as.matrix(e6)),rowSds(as.matrix(e10)))

means <- data.frame(data$time/60,rowMeans(ancestor),rowMeans(control),rowMeans(e6),rowMeans(e10))
names(means) <- c("time", "Ancestor", "C Replicates", "M Replicates","H Replicates")
data2 <- melt(means, id.vars="time")
data2 <- cbind(data2, sd)

names(data2) <- c("Time","Treatment","OD","sd")
data2$Treatment <- factor(data2$Treatment, levels = c("Ancestor", "C Replicates", "M Replicates","H Replicates"))

e <- ggplot(data2, aes(x=Time, y=OD,col=Treatment)) + ylim(-1,.25) +
  geom_point() +  geom_errorbar( aes(ymin = OD-sd, ymax = OD+sd),width = 0.2) + ylab('log(OD)') +  xlab('Time (hr)')  + scale_color_manual(breaks=c("Ancestor", "C Replicates", "M Replicates","H Replicates"), values =c("black","#4daf4a","#ff7f00","#e41a1c"))  + ggtitle("E. Growth Curves in 10% Ethanol YPD") + theme( axis.text=element_text(size=10),
                                                                                                                                                                                                                                                                                                                            axis.title=element_text(size=10,face="bold")) + theme(legend.title = element_text(size=10, face="bold")) + theme(legend.position="bottom") 

#doubling time 
library(growthcurver)
library(growthcurver)
data <- read.table("High_ethanol_2021.txt", header = TRUE)
data <- data[2:nrow(data),]
data$time <- data$time/60
fit_all <- SummarizeGrowthByPlate(data, bg_correct = "blank")
anc <- c(fit_all[1:9,6], NA)
doubling_time <-  data.frame(anc,fit_all[10:19,6],fit_all[20:29,6],fit_all[30:39,6])
Treatment <- c(rep("Ancestor", nrow(doubling_time)),rep("C Replicates", nrow(doubling_time)),rep("M Replicates", nrow(doubling_time)),rep("H Replicates", nrow(doubling_time)))
DD <-   c(anc,fit_all[10:19,6],fit_all[20:29,6],fit_all[30:39,6])
data3 <- data.frame(Treatment, DD)

data3$Treatment <- factor(data3$Treatment, levels = c("Ancestor", "C Replicates", "M Replicates","H Replicates"))
f <- ggplot(data3, aes(x=Treatment, y=DD, fill=Treatment)) +
  ylab('Doubling Time (hr)') +  ggtitle("F. Doubling Time in 10% Ethanol YPD") +
  geom_boxplot()  + scale_fill_manual(breaks=c("Ancestor", "C Replicates", "M Replicates","H Replicates"), values =c("black", "#4daf4a","#ff7f00","#e41a1c")) + ylim(1,4) + theme( axis.text=element_text(size=10),
                                                                                                                                                                                   axis.title=element_text(size=10,face="bold")) + theme(legend.title = element_text(size=15, face="bold")) + theme(legend.position="none")

grid.arrange(a,b,c,d,e,f,ncol=2)


