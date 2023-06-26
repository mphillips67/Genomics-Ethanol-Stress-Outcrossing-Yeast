
library(reshape2)
library("matrixStats")
library(gridExtra)

#plain YPD
#Doubling time comparison
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

tapply(data3$DD,data3$Treatment,mean, na.rm=TRUE)
kruskal.test(DD ~ Treatment, data=data3)
pairwise.wilcox.test(data3$DD, data3$Treatment, p.adjust.method = "BH")

#Carrying capacity comparison
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

tapply(data3$DD,data3$Treatment,mean, na.rm=TRUE)
kruskal.test(DD ~ Treatment, data=data3)
pairwise.wilcox.test(data3$DD, data3$Treatment, p.adjust.method = "BH")


#6% ethanol ypd
#Doubling time comparison
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
kruskal.test(DD ~ Treatment, data=data3)
pairwise.wilcox.test(data3$DD, data3$Treatment, p.adjust.method = "BH")
tapply(data3$DD,data3$Treatment,mean, na.rm=TRUE)

#Carrying capacity comparison
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

kruskal.test(DD ~ Treatment, data=data3)
pairwise.wilcox.test(data3$DD, data3$Treatment, p.adjust.method = "BH")
tapply(data3$DD,data3$Treatment,mean, na.rm=TRUE)


#10% ethaol
#Doubling time comparison
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
kruskal.test(DD ~ Treatment, data=data3)
pairwise.wilcox.test(data3$DD, data3$Treatment, p.adjust.method = "BH")
tapply(data3$DD,data3$Treatment,mean, na.rm=TRUE)


#Carrying capacity comparison
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

kruskal.test(DD ~ Treatment, data=data3)
pairwise.wilcox.test(data3$DD, data3$Treatment, p.adjust.method = "BH")
tapply(data3$DD,data3$Treatment,mean, na.rm=TRUE)




