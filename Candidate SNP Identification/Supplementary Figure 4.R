setwd("~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/Modifed_CMH")

library(gridExtra)
library(ggman)
library(ggpubr)

#make one p-value table and SNPs of interest table
high_ethanol <- read.table('CMH_Results_High_1v15.txt', header = TRUE, col.names = c("CHR","BP","P_H"))
mod_ethanol <- read.table('CMH_Results_Moderate_1v15.txt', header = TRUE, col.names = c("CHR","BP","P_M"))
control <- read.table('CMH_Results_Controls_1v15.txt', header = TRUE, col.names = c("CHR","BP","P_C"))
control$id <- 1:nrow(control)
temp <- merge(high_ethanol, control, by = c("CHR","BP"))
data <- merge(mod_ethanol,temp)
data <- data[order(data$id),]
data$CHR <- as.character(data$CHR)
old <- as.vector(unique(data$CHR))
replacement <- as.numeric(gsub("C","",old))
for (i in 1:16){
  data[data==old[i]] <- replacement[i]
}
SNP <- paste("SNP_", 1:nrow(data),sep="")
data$CHR <- as.numeric(data$CHR)
data <- cbind(SNP,data)

#adjust p-values
data$P_M <- p.adjust(data$P_M, method='fdr')
data$P_H <- p.adjust(data$P_H, method='fdr')
data$P_C <- p.adjust(data$P_C, method='fdr')

#define sig thresholds for all plots for all plots
alpha_threshold <- -log10(0.005)

#High ethanol c1 vs 15
p1 <- ggman(data, snp = "SNP", bp = "BP", chrom = "CHR", ymax = 30, pointSize = .3, pvalue = "P_H", title = "A. High Ethanol Treatment",sigLine = alpha_threshold, relative.positions = TRUE, lineColour = "black",xlabel = "Chromosome", ylabel = "-log10(adj. P-value)")

#Moderate ethanol c1 vs 15
p2 <- ggman(data, snp = "SNP", bp = "BP", chrom = "CHR", ymax = 30, pointSize = .3, pvalue = "P_M", title = "B. Moderate Ethanol Treatment",sigLine = alpha_threshold, relative.positions = TRUE, lineColour = "black",xlabel = "Chromosome", ylabel = "-log10(adj. P-value)")

#control
p3 <- ggman(data, snp = "SNP", bp = "BP", chrom = "CHR", ymax = 30, pointSize = .3, pvalue = "P_C", title = "C. Control Treatment",sigLine = alpha_threshold, relative.positions = TRUE, lineColour = "black",xlabel = "Chromosome", ylabel = "-log10(adj. P-value)")

grid.arrange(p1,p2,p3,nrow=3)

