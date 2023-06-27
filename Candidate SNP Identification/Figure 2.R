setwd("~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/Category Plots")
library(gridExtra)
library(ggman)
library(ggpubr)
#read in sig SNP category and change CHR labels to numbers
outcrossing_snps <- read.table("~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/Outcrossing Snps/Sig_SNPs_Outcrossing.txt", header =TRUE,col.names = c("CHR","BP"))
outcrossing_snps$CHR <- as.character(outcrossing_snps$CHR)
old <- as.vector(unique(outcrossing_snps$CHR))
replacement <- as.numeric(gsub("C","",old))
for (i in 1:16){
  outcrossing_snps[outcrossing_snps==old[i]] <- replacement[i]
}
outcrossing_snps$CHR  <- as.numeric(outcrossing_snps$CHR )

ethanol_shared <- read.table("~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/Ethanol SNPs/Sig_SNPs_Shared_High_Mod_Overlap.txt", header =TRUE,col.names = c("CHR","BP"))
ethanol_shared$CHR <- as.character(ethanol_shared$CHR)
old <- as.vector(unique(ethanol_shared$CHR))
replacement <- as.numeric(gsub("C","",old))
for (i in 1:16){
  ethanol_shared[ethanol_shared==old[i]] <- replacement[i]
}
ethanol_shared$CHR  <- as.numeric(ethanol_shared$CHR )

high_only <- read.table("~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/Ethanol SNPs/High_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
high_only$CHR <- as.character(high_only$CHR)
old <- as.vector(unique(high_only$CHR))
replacement <- as.numeric(gsub("C","",old))
for (i in 1:16){
  high_only[high_only==old[i]] <- replacement[i]
}
high_only$CHR  <- as.numeric(high_only$CHR )

mod_only <- read.table("~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/Ethanol SNPs/Mod_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
mod_only$CHR <- as.character(mod_only$CHR)
old <- as.vector(unique(mod_only$CHR))
replacement <- as.numeric(gsub("C","",old))
for (i in 1:16){
  mod_only[mod_only==old[i]] <- replacement[i]
}
mod_only$CHR  <- as.numeric(mod_only$CHR )

con_only <- read.table("~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/Ethanol SNPs/Control_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
con_only$CHR <- as.character(con_only$CHR)
old <- as.vector(unique(con_only$CHR))
replacement <- as.numeric(gsub("C","",old))
for (i in 1:16){
  con_only[con_only==old[i]] <- replacement[i]
}
con_only$CHR  <- as.numeric(con_only$CHR )

#make one p-value table and SNPs of interest table
high_ethanol <- read.table('~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/CMH/CMH_Results_High_1v15.txt', header = TRUE, col.names = c("CHR","BP","P_H"))
mod_ethanol <- read.table('~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/CMH/CMH_Results_Moderate_1v15.txt', header = TRUE, col.names = c("CHR","BP","P_M"))
control <- read.table('~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/CMH/CMH_Results_Controls_1v15.txt', header = TRUE, col.names = c("CHR","BP","P_C"))
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

outcrossing_snps_highlight <- merge(data, outcrossing_snps, by =c("CHR","BP"))
outcrossing_snps_highlight$group <- rep("General Laboratory Selection",nrow(outcrossing_snps_highlight))

ethanol_snps_highlight <- merge(ethanol_shared,data, by =c("CHR","BP"))
ethanol_snps_highlight$group <- rep("General Ethanol Selection",nrow(ethanol_snps_highlight))

high_snp_highlight <- merge(high_only, data, by =c("CHR","BP"))
high_snp_highlight$group <- rep("High Ethanol Specific",nrow(high_snp_highlight))

mod_snp_highlight <- merge(mod_only, data, by =c("CHR","BP"))
mod_snp_highlight$group <- rep("Moderate Ethanol Specific",nrow(mod_snp_highlight))

con_snp_highlight <- merge(con_only, data, by =c("CHR","BP"))
con_snp_highlight$group <- rep("Control Specific",nrow(con_snp_highlight))

highlight_high <- rbind(ethanol_snps_highlight, outcrossing_snps_highlight, high_snp_highlight)

highlight_mod <- rbind(ethanol_snps_highlight, outcrossing_snps_highlight, mod_snp_highlight)

highlight_con <- rbind(outcrossing_snps_highlight, con_snp_highlight)

highlight <- rbind(ethanol_snps_highlight,outcrossing_snps_highlight, high_snp_highlight, mod_snp_highlight, con_snp_highlight )


#High ethanol c1 vs 15
alpha_high <- -log10(quantile(read.table('~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/CMH/Permutations_CMH/Perms_High_1v15.txt')[,1], .005))
p1 <- ggman(data, snp = "SNP", bp = "BP", chrom = "CHR", ymax = 215, pointSize = .3, pvalue = "P_H", title = "A. High Ethanol Treatment",sigLine = alpha_high, relative.positions = TRUE, lineColour = "black",xlabel = "Chromosome")
high <-ggmanHighlightGroup(p1, highlightDfm = highlight_high, snp = "SNP", group = "group", 
                    size = 1, legend.title = "Groups",legend.remove = FALSE) + scale_fill_manual(breaks=c( "General Laboratory Selection", "General Ethanol Selection", "High Ethanol Specific"), values =c("#377eb8","#984ea3","#e41a1c")) + guides(fill=guide_legend(title="Response Type",override.aes = list(size = 5))) + theme(legend.text = element_text(size = 10),legend.title = element_text(size = 15)) 

#Moderate ethanol c1 vs 15
alpha_mod <- -log10(quantile(read.table('~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/CMH/Permutations_CMH/Perms_Moderate_1v15.txt')[,1], .005))
p2 <- ggman(data, snp = "SNP", bp = "BP", chrom = "CHR", ymax = 215, pointSize = .3, pvalue = "P_M", title = "B. Moderate Ethanol Treatment",sigLine = alpha_mod, relative.positions = TRUE, lineColour = "black",xlabel = "Chromosome")
mod<-ggmanHighlightGroup(p2, highlightDfm = highlight_mod, snp = "SNP", group = "group", 
                    size = 0.5, legend.title = "Groups", legend.remove = FALSE) + scale_fill_manual(breaks=c( "General Laboratory Selection", "General Ethanol Selection","Moderate Ethanol Specific"), values =c("#377eb8","#984ea3","#ff7f00")) + guides(fill=guide_legend(title="Response Type",override.aes = list(size = 5))) + theme(legend.text = element_text(size = 10),legend.title = element_text(size = 15)) 

#control
alpha_cont <- -log10(quantile(read.table('~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/CMH/Permutations_CMH/Perms_Control_1v15.txt')[,1], .005))
p3 <- ggman(data, snp = "SNP", bp = "BP", chrom = "CHR", ymax = 215, pointSize = .3, pvalue = "P_C", title = "C. Control Treatment",sigLine = alpha_cont, relative.positions = TRUE, lineColour = "black",xlabel = "Chromosome")
cont<-ggmanHighlightGroup(p3, highlightDfm = highlight_con, snp = "SNP", group = "group", 
                       size = 0.5, legend.title = "Groups", legend.remove = FALSE) + scale_fill_manual(breaks=c("General Laboratory Selection", "Control Specific"), values =c("#377eb8","#4daf4a")) + guides(fill=guide_legend(title="Response Type",override.aes = list(size = 5)))  + theme(legend.text = element_text(size = 10),legend.title = element_text(size = 15)) 

grid.arrange(high,mod,cont,nrow=3)

