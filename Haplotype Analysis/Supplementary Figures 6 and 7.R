setwd("~/Dropbox/Ethanol expression project/Genomic_analysis/Haplotypes/")
library(gridExtra)
library(ggman)
library(ggpubr)
library(reshape2)
library(stringr)
library(tidyverse)


###Control populations
## process tables with hap data for C1 and C15
control_c1 <-  read.table("Control//Control_wk1_haps.txt", header = TRUE)
founders <- as.vector(str_split_fixed(control_c1[1,5], ";",12))
control_c1_freqs <- cbind(control_c1[,1:3],colsplit(control_c1$adjfounderfreqs, ";", founders))
for (i in 4:ncol(control_c1_freqs)){
  control_c1_freqs <- subset(control_c1_freqs, control_c1_freqs[i] >= 0 & control_c1_freqs[i] <= 1 )
}

control_c15 <-  read.table("Control//Control_wk15_haps.txt", header = TRUE)
founders <- as.vector(str_split_fixed(control_c15[1,5], ";",12))
control_c15_freqs <- cbind(control_c15[,1:3],colsplit(control_c15$adjfounderfreqs, ";", founders))
for (i in 4:ncol(control_c15_freqs)){
  control_c15_freqs <- subset(control_c15_freqs, control_c15_freqs[i] >= 0 & control_c15_freqs[i] <= 1 )
}
control_names_c1 <- unique(control_c1_freqs$poolroot)
control_names_c15 <- unique(control_c15_freqs$poolroot)

haps <- NULL 
#counts
for (i in 1:length(control_names_c1)){
  temp_c1 <- subset(control_c1_freqs, poolroot == control_names_c1[i])
  temp_c15 <- subset(control_c15_freqs, poolroot == control_names_c15[i])
  haps <- rbind(haps, temp_c1[,2:3], temp_c15[,2:3])
}

## calculate D for each replicate 
control_D <- unique(control_c1_freqs[,2:3])
for (i in 1:length(control_names_c1)){
  temp_c1 <- subset(control_c1_freqs, poolroot == control_names_c1[i])
  temp_c15 <- subset(control_c15_freqs, poolroot == control_names_c15[i])
  temp <- merge(temp_c1[,2:15],temp_c15[2:15], by = c("chr","pos"), sort = FALSE)
  temp_D <- (sqrt((rowSums((temp[,3:14] -  temp[,15:26])^2))/12))*100 #calculate %D
  temp_D <- cbind(temp[,1:2], temp_D)
  control_D <- merge(control_D, temp_D, by = c("chr","pos"), sort = FALSE)
  j <- 2+i
  colnames(control_D)[j] <- control_names_c1[i]
}
##calculate average D
control_D$average_D <- rowMeans(control_D[,3:22])
control_D$SD <- apply(control_D[,3:22],1,sd)

##generate plot
SNP <- paste("SNP_", 1:nrow(control_D),sep="")
control_D <- cbind(SNP, control_D)
old <- as.vector(unique(control_D$chr))
replacement <- as.numeric(gsub("C","",old))
for (i in 1:16){
  control_D[control_D==old[i]] <- replacement[i]
}
control_D$chr <- as.numeric(control_D$chr)
p1 <- ggman(control_D, snp = "SNP", bp = "pos", chrom = "chr", ymax = 25, pointSize = .8, pvalue = "average_D", title = "A. Control Treatment",sigLine = NA, relative.positions = TRUE, lineColour = "black",xlabel = "Chromosome",  logTransform = FALSE, ylabel = "Average %D")


###Moderate Populations
moderate_c1 <-  read.table("Moderate/Mod_wk1_haps.txt", header = TRUE)
founders <- as.vector(str_split_fixed(moderate_c1[1,5], ";",12))
moderate_c1_freqs <- cbind(moderate_c1[,1:3],colsplit(moderate_c1$adjfounderfreqs, ";", founders))
for (i in 4:ncol(moderate_c1_freqs)){
  moderate_c1_freqs <- subset(moderate_c1_freqs, moderate_c1_freqs[i] >= 0 & moderate_c1_freqs[i] <= 1 )
}

moderate_c15 <-  read.table("Moderate/Mod_wk15_haps.txt", header = TRUE)
founders <- as.vector(str_split_fixed(moderate_c15[1,5], ";",12))
moderate_c15_freqs <- cbind(moderate_c15[,1:3],colsplit(moderate_c15$adjfounderfreqs, ";", founders))
for (i in 4:ncol(moderate_c15_freqs)){
  moderate_c15_freqs <- subset(moderate_c15_freqs, moderate_c15_freqs[i] >= 0 & moderate_c15_freqs[i] <= 1 )
}

moderate_names_c1 <- unique(moderate_c1_freqs$poolroot)
moderate_names_c15 <- unique(moderate_c15_freqs$poolroot)

#counts
for (i in 1:length(control_names_c1)){
  temp_c1 <- subset(moderate_c1_freqs, poolroot == moderate_names_c1[i])
  temp_c15 <- subset(moderate_c15_freqs, poolroot == moderate_names_c15[i])
  haps <- rbind(haps, temp_c1[,2:3], temp_c15[,2:3])
}

## calculate D for each replicate 
moderate_D <- unique(moderate_c1_freqs[,2:3])
for (i in 1:length(moderate_names_c1)){
  temp_c1 <- subset(moderate_c1_freqs, poolroot == moderate_names_c1[i])
  temp_c15 <- subset(moderate_c15_freqs, poolroot == moderate_names_c15[i])
  temp <- merge(temp_c1[,2:15],temp_c15[2:15], by = c("chr","pos"), sort = FALSE)
  temp_D <- (sqrt((rowSums((temp[,3:14] -  temp[,15:26])^2))/12))*100 #calculate %D
  temp_D <- cbind(temp[,1:2], temp_D)
  moderate_D <- merge(moderate_D, temp_D, by = c("chr","pos"), sort = FALSE)
  j <- 2+i
  colnames(moderate_D)[j] <- moderate_names_c1[i]
}
##calculate average D
moderate_D$average_D <- rowMeans(moderate_D[,3:22])

##generate plot
SNP <- paste("SNP_", 1:nrow(moderate_D),sep="")
moderate_D <- cbind(SNP, moderate_D)
old <- as.vector(unique(moderate_D$chr))
replacement <- as.numeric(gsub("C","",old))
for (i in 1:16){
  moderate_D[moderate_D==old[i]] <- replacement[i]
}
moderate_D$chr <- as.numeric(moderate_D$chr)
p2 <- ggman(moderate_D, snp = "SNP", bp = "pos", chrom = "chr", ymax = 25, pointSize = .8, pvalue = "average_D", title = "B. Moderate Ethanol Treatment",sigLine = NA, relative.positions = TRUE, lineColour = "black",xlabel = "Chromosome",  logTransform = FALSE, ylabel = "Average %D")


###High populations
## process tables with hap data for C1 and C15
high_c1 <-  read.table("High/High_wk1_haps.txt", header = TRUE)
founders <- as.vector(str_split_fixed(high_c1[1,5], ";",12))
high_c1_freqs <- cbind(high_c1[,1:3],colsplit(high_c1$adjfounderfreqs, ";", founders))
for (i in 4:ncol(high_c1_freqs)){
  high_c1_freqs <- subset(high_c1_freqs, high_c1_freqs[i] >= 0 & high_c1_freqs[i] <= 1 )
}

high_c15 <-  read.table("High/High_wk15_haps.txt", header = TRUE)
founders <- as.vector(str_split_fixed(high_c15[1,5], ";",12))
high_c15_freqs <- cbind(high_c15[,1:3],colsplit(high_c15$adjfounderfreqs, ";", founders))
for (i in 4:ncol(high_c15_freqs)){
  high_c15_freqs <- subset(high_c15_freqs, high_c15_freqs[i] >= 0 & high_c15_freqs[i] <= 1 )
}
high_names_c1 <- unique(high_c1_freqs$poolroot)
high_names_c15 <- unique(high_c15_freqs$poolroot)

#counts
for (i in 1:length(control_names_c1)){
  temp_c1 <- subset(high_c1_freqs, poolroot == high_names_c1[i])
  temp_c15 <- subset(high_c15_freqs, poolroot == high_names_c15[i])
  haps <- rbind(haps, temp_c1[,2:3], temp_c15[,2:3])
}


## calculate D for each replicate 
high_D <- unique(high_c1_freqs[,2:3])
for (i in 1:length(high_names_c1)){
  temp_c1 <- subset(high_c1_freqs, poolroot == high_names_c1[i])
  temp_c15 <- subset(high_c15_freqs, poolroot == high_names_c15[i])
  temp <- merge(temp_c1[,2:15],temp_c15[2:15], by = c("chr","pos"), sort = FALSE)
  temp_D <- (sqrt((rowSums((temp[,3:14] -  temp[,15:26])^2))/12))*100 #calculate %D
  temp_D <- cbind(temp[,1:2], temp_D)
  high_D <- merge(high_D, temp_D, by = c("chr","pos"), sort = FALSE)
  j <- 2+i
  colnames(high_D)[j] <- high_names_c1[i]
}
##calculate average D
high_D$average_D <- rowMeans(high_D[,3:22])

##generate plot
SNP <- paste("SNP_", 1:nrow(high_D),sep="")
high_D <- cbind(SNP, high_D)
old <- as.vector(unique(high_D$chr))
replacement <- as.numeric(gsub("C","",old))
for (i in 1:16){
  high_D[high_D==old[i]] <- replacement[i]
}
high_D$chr <- as.numeric(high_D$chr)
p3 <- ggman(high_D, snp = "SNP", bp = "pos", chrom = "chr", ymax = 25, pointSize = .8, pvalue = "average_D", title = "C. High Ethanol Treatment",sigLine = NA, relative.positions = TRUE, lineColour = "black",xlabel = "Chromosome",  logTransform = FALSE, ylabel = "Average %D")


# plot figs together
grid.arrange(p1,p2,p3,nrow=3)


#make desnity plots of average D
#control
(sum(control_D$average_D <= 10) / length(control_D$average_D)) * 100 #10 or less
(sum(control_D$average_D >= 20) / length(control_D$average_D)) * 100 #15 or more

a <- ggplot(control_D, aes(x = average_D)) +
  geom_density(fill = "blue", color = "blue",alpha = 0.7) +
  labs(title = "A. Control Treatment",
       x = "Average %D",
       y = "Density") +
  theme_minimal() +
  ylim(0, .25)+
  xlim(0, 21)



#moderate
(sum(moderate_D$average_D <= 10) / length(moderate_D$average_D)) * 100 #10 or less
(sum(moderate_D$average_D >= 20) / length(moderate_D$average_D)) * 100 #15 or more

b <- ggplot(moderate_D, aes(x = average_D)) +
  geom_density(fill = "blue", color = "blue",alpha = 0.7) +
  labs(title = "B. Moderate Ethanol Treatment",
       x = "Average %D",
       y = "Density") +
  theme_minimal()+
  ylim(0, .25)+
  xlim(0, 21)



#high
(sum(high_D$average_D <= 10) / length(high_D$average_D)) * 100 #10 or less
(sum(high_D$average_D >= 20) / length(high_D$average_D)) * 100 #15 or more

c <- ggplot(high_D, aes(x = average_D)) +
  geom_density(fill = "blue", color = "blue",alpha = 0.7) +
  labs(title = "C. High Ethanol Treatment",
       x = "Average %D",
       y = "Density") +
  theme_minimal()+
  ylim(0, .25)+
  xlim(0, 21)


grid.arrange(a,b,c,nrow=3)





qqnorm(control_D$average_D)
qqline(control_D$average_D)

qqnorm(moderate_D$average_D)
qqline(moderate_D$average_D)

qqnorm(high_D$average_D)
qqline(high_D$average_D)
