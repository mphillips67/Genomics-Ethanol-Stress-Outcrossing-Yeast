
#define sig thresholds for all plots for all lits
alpha_threshold <- 0.005

#control 1 vs 15
Control_1v15 <- read.table('CMH_Results_Controls_1v15.txt', header = TRUE)
Control_1v15$pval <- p.adjust(Control_1v15$pval, method = "fdr")
Sig_Control_1v15 <- subset(Control_1v15, pval <= alpha_threshold)
write.table(Sig_Control_1v15,"Sig_SNPs_Control_1v15.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

#moderate 1 vs 15
Moderate_1v15 <- read.table('CMH_Results_Moderate_1v15.txt', header = TRUE)
Moderate_1v15$pval <- p.adjust(Moderate_1v15$pval, method = "fdr")
Sig_Moderate_1v15 <- subset(Moderate_1v15, pval <= alpha_threshold)
write.table(Sig_Moderate_1v15,"Sig_SNPs_Moderate_1v15.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

#high 1 vs 15
High_1v15 <- read.table('CMH_Results_High_1v15.txt', header = TRUE)
High_1v15$pval <- p.adjust(High_1v15$pval, method = "fdr")
Sig_High_1v15 <- subset(High_1v15, pval <= alpha_threshold)
write.table(Sig_High_1v15,"Sig_SNPs_High_1v15.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")
