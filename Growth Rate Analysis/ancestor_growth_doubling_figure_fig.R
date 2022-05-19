
library(ggplot2)
library(reshape2)
library("matrixStats")

#plain ypd
data <- read.table("Eth_48_hr_timepoint_0.txt", header = TRUE)

#adjust for different in media color (autoclaved and vacuum filtered)
data[,3:22] <- data[,3:22] - .04
data[,23:42] <- data[,23:42] - .015

data[,2:ncol(data)] <- log10(data[,2:ncol(data)])
#data <- data[1:62,]
time <- data$time/60
control <- data[,3:22] 
e6 <- data[,23:42]
e10 <- data[,43:62]

control_mean <- cbind(time, rowMeans(control), rep("Control", nrow(control)))
e6_mean <- cbind(time, rowMeans(e6), rep("6% ethanol", nrow(e6)))
e10_mean <- cbind(time, rowMeans(e10), rep("10% ethanol", nrow(e10)))
sd <- c(rowSds(as.matrix(control)),rowSds(as.matrix(e6)),rowSds(as.matrix(e10)))
means <- data.frame(data$time/60,rowMeans(control),rowMeans(e6),rowMeans(e10))
names(means) <- c("time", "YPD", "6% Ethanol YPD","10% Ethanol YPD")
data2 <- melt(means, id.vars="time")
data2 <- cbind(data2, sd)

names(data2) <- c("Time","Media","OD","sd")

p <- ggplot(data2, aes(x=Time, y=OD,col=Media)) +  ylim(-1,.3) +
  geom_point() +  geom_errorbar( aes(ymin = OD-sd, ymax = OD+sd),width = 0.2) + ylab('log(OD)') +  xlab('Time (hr)')  + scale_color_manual(breaks=c("YPD", "6% Ethanol YPD", "10% Ethanol YPD"), values =c("#4daf4a","#ff7f00","#e41a1c"))  + ggtitle("A. Growth Curves Ancestral Population") + theme( axis.text=element_text(size=10),
                                                                                                                                                                                                                                                                                                        axis.title=element_text(size=10,face="bold")) + theme(legend.title = element_text(size=10, face="bold")) 

#doubling time 
library(growthcurver)
data <- read.table("Eth_48_hr_timepoint_0.txt", header = TRUE)
data$time <- data$time/60
fit_all <- SummarizeGrowthByPlate(data, bg_correct = "blank")

doubling_time <-  data.frame(fit_all[1:20,6],fit_all[21:40,6],fit_all[41:60,6])
Media <- c(rep("YPD", nrow(doubling_time)),rep("6% Ethanol YPD", nrow(doubling_time)),rep("10% Ethanol YPD", nrow(doubling_time)))
DD <-  c(fit_all[1:20,6],fit_all[21:40,6],fit_all[41:60,6])

data3 <- data.frame(Media, DD)

data3$Media <- factor(data3$Media, levels = c("YPD", "6% Ethanol YPD", "10% Ethanol YPD"))
q <- ggplot(data3, aes(x=Media, y=DD, fill=Media)) +
  ylab('Doubling Time (hr)') +  ggtitle("B. Doubling Time Ancestral Population") +
  geom_boxplot()  + scale_fill_manual(breaks=c("YPD", "6% Ethanol YPD", "10% Ethanol YPD"), values =c("#4daf4a","#ff7f00","#e41a1c")) + ylim(1,3) + theme( axis.text=element_text(size=10),
                                                                                                                                                           axis.title=element_text(size=10,face="bold")) + theme(legend.title = element_text(size=15, face="bold")) 


library(gridExtra)
  grid.arrange(p,q,ncol=2)


