
library(ggplot2)
library(ggpubr)
library(gginference)
#SNP table
data <- read.table("~/Dropbox/Ethanol expression project/Genomic_analysis/Prep_quality_checks/SNPtable_raw.txt", header = TRUE)

#All SNP categories
#read in category lists
high_only <- read.table("High_only_sig_snps.txt", header =TRUE,col.names = c("chr","pos"))
mod_only <- read.table("Mod_only_sig_snps.txt", header =TRUE,col.names = c("chr","pos"))
con_only <- read.table("Control_only_sig_snps.txt", header =TRUE,col.names = c("chr","pos"))
general_lab  <- read.table("Sig_SNPs_Outcrossing.txt", header =TRUE,col.names = c("chr","pos"))
general_eth <- read.table("Sig_SNPs_Shared_High_Mod_Overlap.txt", header =TRUE,col.names = c("chr","pos"))

#combined SNP categories for ethanol 
ethnaol_cands <- rbind(high_only, general_eth)


#treatment + tps
control <-c(paste("0",1:9,sep=""),10,12:21)
moderate <- c(25:30,32:33,35:38,41:48)
high <- c(49:56,59:61,64:72)

#lgeneral ab candidates plot
#chang in high 
Sig_SNPS_High <- merge(data, general_lab, by = c("chr","pos"))
maf_week1 <- paste("alt_ETH_rep", high,"_wk01",sep="")
maf_week1_table <- Sig_SNPS_High[,c(maf_week1)]
maf_week15 <- paste("alt_ETH_rep", high,"_wk15",sep="")
maf_week15_table <- Sig_SNPS_High[,c(maf_week15)]
change_maf_high <- abs(maf_week15_table-maf_week1_table)
mean_change_high <- colMeans(change_maf_high)

#chang in mod 
Sig_SNPS_Mod <- merge(data, general_lab, by = c("chr","pos"))
maf_week1 <- paste("alt_ETH_rep", moderate,"_wk01",sep="")
maf_week1_table <- Sig_SNPS_Mod[,c(maf_week1)]
maf_week15 <- paste("alt_ETH_rep", moderate,"_wk15",sep="")
maf_week15_table <- Sig_SNPS_Mod[,c(maf_week15)]
change_maf_moderate<- abs(maf_week15_table-maf_week1_table)
mean_change_moderate <- colMeans(change_maf_moderate)

#change in control
Sig_SNPS_Control <- merge(data, general_lab, by = c("chr","pos"))
maf_week1 <- paste("alt_ETH_rep", control,"_wk01",sep="")
maf_week1_table <- Sig_SNPS_Control[,c(maf_week1)]
maf_week15 <- paste("alt_ETH_rep", control,"_wk15",sep="")
maf_week15_table <- Sig_SNPS_Control[,c(maf_week15)]
change_maf_control <- abs(maf_week15_table-maf_week1_table)
mean_change_control <- colMeans(change_maf_control)

#boxplot 
delta_maf <- c(mean_change_control, mean_change_moderate,mean_change_high)
Treatment <- c(rep("Control", length(mean_change_control)),rep("Moderate Ethanol", length(mean_change_moderate)),rep("High Ethanol", length(mean_change_high)))
data_frame <- data.frame(Treatment, delta_maf)
data_frame$Treatment <- factor(data_frame$Treatment, levels = c("Control", "Moderate Ethanol", "High Ethanol"))

my_comparisons <- list(c("High Ethanol", "Moderate Ethanol"), c("Control", "High Ethanol"), c("Control", "Moderate Ethanol"))

a <- ggplot(data_frame, aes(x=Treatment, y=delta_maf, fill=Treatment)) +
  ylab('Change in SNP Frequency') + 
  geom_boxplot() +  stat_compare_means(comparisons = my_comparisons, method="wilcox.test", label = "p.signif", size=3, vjust=0.9)  + scale_fill_manual(breaks=c("Control", "Moderate Ethanol", "High Ethanol"), values =c("#4daf4a","#ff7f00","#e41a1c")) +  stat_compare_means(label.y = 0.05, label.x =2, size=5) +  ylim(0,.5) + ggtitle("A.  General Lab Selection Candidates") + theme_grey(base_size = 12) + theme(legend.position="none")


#ethanol candidate plots 
#chang in high 
Sig_SNPS_High <- merge(data, ethnaol_cands, by = c("chr","pos"))
maf_week1 <- paste("alt_ETH_rep", high,"_wk01",sep="")
maf_week1_table <- Sig_SNPS_High[,c(maf_week1)]
maf_week15 <- paste("alt_ETH_rep", high,"_wk15",sep="")
maf_week15_table <- Sig_SNPS_High[,c(maf_week15)]
change_maf_high <- abs(maf_week15_table-maf_week1_table)
mean_change_high <- colMeans(change_maf_high)

#chang in mod 
Sig_SNPS_Mod <- merge(data, ethnaol_cands, by = c("chr","pos"))
maf_week1 <- paste("alt_ETH_rep", moderate,"_wk01",sep="")
maf_week1_table <- Sig_SNPS_Mod[,c(maf_week1)]
maf_week15 <- paste("alt_ETH_rep", moderate,"_wk15",sep="")
maf_week15_table <- Sig_SNPS_Mod[,c(maf_week15)]
change_maf_moderate<- abs(maf_week15_table-maf_week1_table)
mean_change_moderate <- colMeans(change_maf_moderate)

#change in control
Sig_SNPS_Control <- merge(data, ethnaol_cands, by = c("chr","pos"))
maf_week1 <- paste("alt_ETH_rep", control,"_wk01",sep="")
maf_week1_table <- Sig_SNPS_Control[,c(maf_week1)]
maf_week15 <- paste("alt_ETH_rep", control,"_wk15",sep="")
maf_week15_table <- Sig_SNPS_Control[,c(maf_week15)]
change_maf_control <- abs(maf_week15_table-maf_week1_table)
mean_change_control <- colMeans(change_maf_control)

#boxplot 
delta_maf <- c(mean_change_control, mean_change_moderate,mean_change_high)
Treatment <- c(rep("Control", length(mean_change_control)),rep("Moderate Ethanol", length(mean_change_moderate)),rep("High Ethanol", length(mean_change_high)))
data_frame <- data.frame(Treatment, delta_maf)
data_frame$Treatment <- factor(data_frame$Treatment, levels = c("Control", "Moderate Ethanol", "High Ethanol"))

my_comparisons <- list(c("High Ethanol", "Moderate Ethanol"), c("Control", "High Ethanol"), c("Control", "Moderate Ethanol"))

c <- ggplot(data_frame, aes(x=Treatment, y=delta_maf, fill=Treatment)) +
  ylab('Change in SNP Frequency') + 
  geom_boxplot()  +  stat_compare_means(comparisons = my_comparisons, method="wilcox.test", label = "p.signif", size=3, vjust=0.9)  + scale_fill_manual(breaks=c("Control", "Moderate Ethanol", "High Ethanol"), values =c("#4daf4a","#ff7f00","#e41a1c")) +  stat_compare_means(label.y = 0.05, label.x =2, size=5) +  ylim(0,.5) + ggtitle("C. General Ethanol Selection Candidates") + theme_grey(base_size = 12) + theme(legend.position="none")


#moderate specific candidates
#chang in high 
Sig_SNPS_High <- merge(data, mod_only, by = c("chr","pos"))
maf_week1 <- paste("alt_ETH_rep", high,"_wk01",sep="")
maf_week1_table <- Sig_SNPS_High[,c(maf_week1)]
maf_week15 <- paste("alt_ETH_rep", high,"_wk15",sep="")
maf_week15_table <- Sig_SNPS_High[,c(maf_week15)]
change_maf_high <- abs(maf_week15_table-maf_week1_table)
mean_change_high <- colMeans(change_maf_high)

#chang in mod 
Sig_SNPS_Mod <- merge(data, mod_only, by = c("chr","pos"))
maf_week1 <- paste("alt_ETH_rep", moderate,"_wk01",sep="")
maf_week1_table <- Sig_SNPS_Mod[,c(maf_week1)]
maf_week15 <- paste("alt_ETH_rep", moderate,"_wk15",sep="")
maf_week15_table <- Sig_SNPS_Mod[,c(maf_week15)]
change_maf_moderate<- abs(maf_week15_table-maf_week1_table)
mean_change_moderate <- colMeans(change_maf_moderate)

#change in control
Sig_SNPS_Control <- merge(data, mod_only, by = c("chr","pos"))
maf_week1 <- paste("alt_ETH_rep", control,"_wk01",sep="")
maf_week1_table <- Sig_SNPS_Control[,c(maf_week1)]
maf_week15 <- paste("alt_ETH_rep", control,"_wk15",sep="")
maf_week15_table <- Sig_SNPS_Control[,c(maf_week15)]
change_maf_control <- abs(maf_week15_table-maf_week1_table)
mean_change_control <- colMeans(change_maf_control)

#boxplot 
delta_maf <- c(mean_change_control, mean_change_moderate,mean_change_high)
Treatment <- c(rep("Control", length(mean_change_control)),rep("Moderate Ethanol", length(mean_change_moderate)),rep("High Ethanol", length(mean_change_high)))
data_frame <- data.frame(Treatment, delta_maf)
data_frame$Treatment <- factor(data_frame$Treatment, levels = c("Control", "Moderate Ethanol", "High Ethanol"))

my_comparisons <- list(c("High Ethanol", "Moderate Ethanol"), c("Control", "High Ethanol"), c("Control", "Moderate Ethanol"))

d <- ggplot(data_frame, aes(x=Treatment, y=delta_maf, fill=Treatment)) +
  ylab('Change in SNP Frequency') + 
  geom_boxplot()  +  stat_compare_means(comparisons = my_comparisons, method="wilcox.test", label = "p.signif",size=3, vjust=1.1)  + scale_fill_manual(breaks=c("Control", "Moderate Ethanol", "High Ethanol"), values =c("#4daf4a","#ff7f00","#e41a1c")) + stat_compare_means(label.y = 0.05, label.x =2, size=5) +  ylim(0,.5) + ggtitle("D. Moderate Ethanol Specific Candidates") + theme_grey(base_size = 12) + theme(legend.position="none")


#control specific plots
#chang in high 
Sig_SNPS_High <- merge(data, con_only, by = c("chr","pos"))
maf_week1 <- paste("alt_ETH_rep", high,"_wk01",sep="")
maf_week1_table <- Sig_SNPS_High[,c(maf_week1)]
maf_week15 <- paste("alt_ETH_rep", high,"_wk15",sep="")
maf_week15_table <- Sig_SNPS_High[,c(maf_week15)]
change_maf_high <- abs(maf_week15_table-maf_week1_table)
mean_change_high <- colMeans(change_maf_high)

#chang in mod 
Sig_SNPS_Mod <- merge(data, con_only, by = c("chr","pos"))
maf_week1 <- paste("alt_ETH_rep", moderate,"_wk01",sep="")
maf_week1_table <- Sig_SNPS_Mod[,c(maf_week1)]
maf_week15 <- paste("alt_ETH_rep", moderate,"_wk15",sep="")
maf_week15_table <- Sig_SNPS_Mod[,c(maf_week15)]
change_maf_moderate<- abs(maf_week15_table-maf_week1_table)
mean_change_moderate <- colMeans(change_maf_moderate)

#change in control
Sig_SNPS_Control <- merge(data, con_only, by = c("chr","pos"))
maf_week1 <- paste("alt_ETH_rep", control,"_wk01",sep="")
maf_week1_table <- Sig_SNPS_Control[,c(maf_week1)]
maf_week15 <- paste("alt_ETH_rep", control,"_wk15",sep="")
maf_week15_table <- Sig_SNPS_Control[,c(maf_week15)]
change_maf_control <- abs(maf_week15_table-maf_week1_table)
mean_change_control <- colMeans(change_maf_control)

#boxplot 
delta_maf <- c(mean_change_control, mean_change_moderate,mean_change_high)
Treatment <- c(rep("Control", length(mean_change_control)),rep("Moderate Ethanol", length(mean_change_moderate)),rep("High Ethanol", length(mean_change_high)))
data_frame <- data.frame(Treatment, delta_maf)
data_frame$Treatment <- factor(data_frame$Treatment, levels = c("Control", "Moderate Ethanol", "High Ethanol"))

my_comparisons <- list(c("High Ethanol", "Moderate Ethanol"), c("Control", "High Ethanol"), c("Control", "Moderate Ethanol"))

b <- ggplot(data_frame, aes(x=Treatment, y=delta_maf, fill=Treatment)) +
  ylab('Change in SNP Frequency') + 
  geom_boxplot()  +  stat_compare_means(comparisons = my_comparisons, method="wilcox.test", label = "p.signif",size=3, vjust=1.1)  + scale_fill_manual(breaks=c("Control", "Moderate Ethanol", "High Ethanol"), values =c("#4daf4a","#ff7f00","#e41a1c")) + stat_compare_means(label.y = 0.05, label.x =2, size=5) +  ylim(0,.5) + ggtitle("B. Control Treatment Specific Candidates") + theme_grey(base_size = 12) + theme(legend.position="none")


library(gridExtra)
grid.arrange(a,b,c,d,ncol=2)

