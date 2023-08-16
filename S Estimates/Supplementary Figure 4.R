
library(ggman)

#read in results and prep for plotting 
data <- read.table("BaitER_results_all.txt", header = TRUE, col.names = c("CHR","BP","Ref","Control","Moderate","High"))
data$CHR <- as.character(data$CHR)
old <- as.vector(unique(data$CHR))
replacement <- as.numeric(gsub("C","",old))
for (i in 1:16){
  data[data==old[i]] <- replacement[i]
}
SNP <- paste("SNP_", 1:nrow(data),sep="")
data$CHR <- as.numeric(data$CHR)
data$SNP <-SNP



#High ethanol 
p1 <- ggman(data, snp = "SNP", bp = "BP", chrom = "CHR", logTransform = FALSE, ymax = 0.2, ymin =-0.02, pointSize = .5, pvalue = "High", title = "A. High Ethanol Treatment",sigLine = NA, ylabel = "Selection Coefficient", relative.positions = TRUE, lineColour = "black",xlabel = "Chromosome")

#Moderate ethanol
p2 <- ggman(data, snp = "SNP", bp = "BP", chrom = "CHR", logTransform = FALSE, ymax = 0.2, ymin = -0.02, pointSize = .5, pvalue = "Moderate", title = "B. Moderate Ethanol Treatment",sigLine = NA, ylabel = "Selection Coefficient", relative.positions = TRUE, lineColour = "black",xlabel = "Chromosome")

#control
p3 <- ggman(data, snp = "SNP", bp = "BP", chrom = "CHR", logTransform = FALSE, ymax = 0.2, ymin = -0.02, pointSize = .5, pvalue = "Control", title = "C. Control Treatment",sigLine = NA, ylabel = "Selection Coefficient", relative.positions = TRUE, lineColour = "black",xlabel = "Chromosome")


grid.arrange(p1,p2,p3,nrow=3)

