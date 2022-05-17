
library(ggplot2)
library(reshape2)
library("matrixStats")

#carrying capacity plain

library(growthcurver)
data <- read.table("Plain_ypd_2021.txt", header = TRUE)
data <- data[2:nrow(data),]
data[,2:ncol(data)] <- log10(data[,2:ncol(data)])
data$time <- data$time/60
fit_all <- SummarizeGrowthByPlate(data, bg_correct = "blank")
anc <- c(fit_all[1:9,2], NA)
doubling_time <-  data.frame(anc,fit_all[10:19,2],fit_all[20:29,2],fit_all[30:39,2])
Treatment <- c(rep("Ancestor", nrow(doubling_time)),rep("C Replicates", nrow(doubling_time)),rep("M Replicates", nrow(doubling_time)),rep("H Replicates", nrow(doubling_time)))
DD <-   c(anc,fit_all[10:19,2],fit_all[20:29,2],fit_all[30:39,2])
data3 <- data.frame(Treatment, DD)
data3$Treatment <- factor(data3$Treatment, levels = c("Ancestor", "C Replicates", "M Replicates","H Replicates"))
a <- ggplot(data3, aes(x=Treatment, y=DD, fill=Treatment)) +
  ylab('Carrying Capacity (log(OD))') +  ggtitle("A. Carrying Capacity in Plain YPD") +
  geom_boxplot()  + scale_fill_manual(breaks=c("Ancestor", "C Replicates", "M Replicates","H Replicates"), values =c("black", "#4daf4a","#ff7f00","#e41a1c")) + ylim(0,2) + theme( axis.text=element_text(size=10),                                                                                                                                                   axis.title=element_text(size=10,face="bold")) + theme(legend.title = element_text(size=15, face="bold"))



#carrying cpacity mod
library(growthcurver)
data <- read.table("Moderate_ethanol_2021.txt", header = TRUE)
data <- data[2:nrow(data),]
data[,2:ncol(data)] <- log10(data[,2:ncol(data)])
data$time <- data$time/60
fit_all <- SummarizeGrowthByPlate(data, bg_correct = "blank")
anc <- c(fit_all[1:9,2], NA)
doubling_time <-  data.frame(anc,fit_all[10:19,2],fit_all[20:29,2],fit_all[30:39,2])
Treatment <- c(rep("Ancestor", nrow(doubling_time)),rep("C Replicates", nrow(doubling_time)),rep("M Replicates", nrow(doubling_time)),rep("H Replicates", nrow(doubling_time)))
DD <-   c(anc,fit_all[10:19,2],fit_all[20:29,2],fit_all[30:39,2])
data3 <- data.frame(Treatment, DD)


data3$Treatment <- factor(data3$Treatment, levels = c("Ancestor", "C Replicates", "M Replicates","H Replicates"))
b <- ggplot(data3, aes(x=Treatment, y=DD, fill=Treatment)) +
  ylab('Carrying Capacity (log(OD))') +  ggtitle("B. Carrying Capacity in 6% Ethanol YPD") +
  geom_boxplot()  + scale_fill_manual(breaks=c("Ancestor", "C Replicates", "M Replicates","H Replicates"), values =c("black", "#4daf4a","#ff7f00","#e41a1c")) + ylim(0,2) + theme( axis.text=element_text(size=10),
                                                                                                                                                                                   axis.title=element_text(size=10,face="bold")) + theme(legend.title = element_text(size=15, face="bold"))





#carrying capacity high
library(growthcurver)
data <- read.table("High_ethanol_2021.txt", header = TRUE)
data <- data[2:nrow(data),]
data[,2:ncol(data)] <- log10(data[,2:ncol(data)])
data$time <- data$time/60
fit_all <- SummarizeGrowthByPlate(data, bg_correct = "blank")
anc <- c(fit_all[1:9,2], NA)
doubling_time <-  data.frame(anc,fit_all[10:19,2],fit_all[20:29,2],fit_all[30:39,2])
Treatment <- c(rep("Ancestor", nrow(doubling_time)),rep("C Replicates", nrow(doubling_time)),rep("M Replicates", nrow(doubling_time)),rep("H Replicates", nrow(doubling_time)))
DD <-   c(anc,fit_all[10:19,2],fit_all[20:29,2],fit_all[30:39,2])
data3 <- data.frame(Treatment, DD)


data3$Treatment <- factor(data3$Treatment, levels = c("Ancestor", "C Replicates", "M Replicates","H Replicates"))
c <- ggplot(data3, aes(x=Treatment, y=DD, fill=Treatment)) +
  ylab('Carrying Capacity (log(OD))') +  ggtitle("C. Carrying Capacity in 10% Ethanol YPD") +
  geom_boxplot()  + scale_fill_manual(breaks=c("Ancestor", "C Replicates", "M Replicates","H Replicates"), values =c("black", "#4daf4a","#ff7f00","#e41a1c")) + ylim(0,2) + theme( axis.text=element_text(size=10),
                                                                                                                                                                                   axis.title=element_text(size=10,face="bold")) + theme(legend.title = element_text(size=15, face="bold"))                                                                                                                                                                            
grid.arrange(a,b,c,ncol=3)




