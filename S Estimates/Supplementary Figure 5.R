library(cowplot)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(grid)
library(ggcorrplot)
library(dplyr)
library(GGally)

#read in results and make S table
control_s <- read.table("Control_BaitER_results.txt", header = TRUE)
mod_s <- read.table("Moderate_BaitER_results.txt", header = TRUE)
high_s <- read.table("High_BaitER_results.txt", header = TRUE)
data <- data.frame(control_s[,1:4], mod_s$sigma, high_s$sigma)
names(data) <- c("CHR", "BP", "REF", "C-populations","M-populations", "H-populations")


#read in category lists
high_only <- read.table("High_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
mod_only <- read.table("Mod_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
control_only <- read.table("Control_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
general_lab  <- read.table("Sig_SNPs_Outcrossing.txt", header =TRUE,col.names = c("CHR","BP"))
general_eth <- read.table("Sig_SNPs_Shared_High_Mod_Overlap.txt", header =TRUE,col.names = c("CHR","BP"))



#genome wide s r2 plot
corr_all <- cor(data[,4:6], use = "complete.obs", method = "pearson") #cor matrix
all <- data[,4:6]
a <- ggpairs(all, title = "A. All Polymorphic Sites")


#general lab 
lab <- na.omit(merge(data, general_lab))[,4:6]
b <- ggpairs(lab, title = "B. General Lab Selection Candidates")

#general ethanol
ethanol <- na.omit(merge(data, general_eth))[,4:6]
c <- ggpairs(ethanol, title = "C. General Ethanol Selection Candidates")

#high ethanol specific 
high <- na.omit(merge(data, high_only))[,4:6]
d <- ggpairs(high, title = "D. High Ethanol Specific Candidates")

#moderate specific
moderate <- na.omit(merge(data, mod_only))[,4:6]
e <- ggpairs(moderate, title = "E. Moderate Ethanol Specific Candidates")

#control specific
control <- na.omit(merge(data, control_only))[,4:6]
f <- ggpairs(control, title = "F. Control Specific Candidates")


plot_grid(
  ggmatrix_gtable(a),
  ggmatrix_gtable(b),
  ggmatrix_gtable(c),
  ggmatrix_gtable(d),
  ggmatrix_gtable(e),
  ggmatrix_gtable(f),
  nrow = 2)
