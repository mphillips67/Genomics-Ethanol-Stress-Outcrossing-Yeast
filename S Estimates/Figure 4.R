#requires category plots - see candidate SNP idenetification directory for how to generate these

library(gridExtra)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(grid)
library(ggcorrplot)

#read in results and make S table
control <- read.table("Control_BaitER_results.txt", header = TRUE)
mod <- read.table("Moderate_BaitER_results.txt", header = TRUE)
high <- read.table("High_BaitER_results.txt", header = TRUE)
data <- data.frame(control[,1:4], mod$sigma, high$sigma)
names(data) <- c("CHR", "BP", "REF", "Control","Moderate", "High")

#read in category lists
high_only <- read.table("High_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
mod_only <- read.table("Mod_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
con_only <- read.table("Control_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))

#genome wide s r2 plot
corr_all <- cor(data[,4:6], use = "complete.obs", method = "pearson")
a <- ggcorrplot(corr_all^2, method = "square", outline.col = "white", lab = TRUE, type = "lower",
                title = expression("A. Genome-wide " ~ italic(s))) +
  scale_fill_gradient(low = "white", high = "#E46726", limits = c(0, 1)) +
  labs(fill = expression(R^2)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

#high only s r2 plot 
high <- na.omit(merge(data, high_only))
corr_high <- cor(high[, 4:6], use = "complete.obs", method = "pearson")
b <- ggcorrplot(corr_high^2, method = "square", outline.col = "white", lab = TRUE, type = "lower",
                title = expression("B. High Ethanol Specific " ~ italic(s))) +
  scale_fill_gradient(low = "white", high = "#E46726", limits = c(0, 1)) +
  labs(fill = expression(R^2)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

#mod only s r2 plot
moderate <- na.omit(merge(data, mod_only))
corr_mod <- cor(moderate[, 4:6], use = "complete.obs", method = "pearson")
c <- ggcorrplot(corr_mod^2, method = "square", outline.col = "white", lab = TRUE, type = "lower",
                title = expression("C. Moderate Ethanol Specific " ~ italic(s))) +
  scale_fill_gradient(low = "white", high = "#E46726", limits = c(0, 1)) +
  labs(fill = expression(R^2)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

#control only s r2 plot
control <- na.omit(merge(data, con_only))
corr_con <- cor(control[, 4:6], use = "complete.obs", method = "pearson")
d <- ggcorrplot(corr_con^2, method = "square", outline.col = "white", lab = TRUE, type = "lower",
                title = expression("D. Control Specific " ~ italic(s))) +
  scale_fill_gradient(low = "white", high = "#E46726", limits = c(0, 1)) +
  labs(fill = expression(R^2)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))


grid.arrange(a,b,c,d,ncol=2)


