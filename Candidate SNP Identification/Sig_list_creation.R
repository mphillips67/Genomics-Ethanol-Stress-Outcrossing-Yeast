#control 1 vs 15
Control_1v15 <- read.table('CMH_Results_Controls_1v15.txt', header = TRUE)
alpha1 <- quantile(read.table('Permutations_CMH/Perms_Control_1v15.txt')[,1],.005)
Sig_Control_1v15 <- subset(Control_1v15, pval <= alpha1)
write.table(Sig_Control_1v15,"Sig_SNPs_Control_1v15.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

#moderate 1 vs 15
Moderate_1v15 <- read.table('CMH_Results_Moderate_1v15.txt', header = TRUE)
alpha2 <- quantile(read.table('Permutations_CMH/Perms_Moderate_1v15.txt')[,1],.005)
Sig_Moderate_1v15 <- subset(Moderate_1v15, pval <= alpha2)
write.table(Sig_Moderate_1v15,"Sig_SNPs_Moderate_1v15.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

#high 1 vs 15
High_1v15 <- read.table('CMH_Results_High_1v15.txt', header = TRUE)
alpha3 <- quantile(read.table('Permutations_CMH/Perms_High_1v15.txt')[,1],.005)
Sig_High_1v15 <- subset(High_1v15, pval <= alpha3)
write.table(Sig_High_1v15,"Sig_SNPs_High_1v15.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

#glmm high vs mod
HvM_GLMM <- read.table('GLMM_Results_HighvMod.txt', header =TRUE)
alpha7 <- quantile(read.table('Permutations_GLMM/Perms_GLMM_HighvMod.txt')[,1], .005)
Sig_HvM_GLMM <- subset(HvM_GLMM, pval <= alpha7)
write.table(Sig_HvM_GLMM,"Sig_SNPs_GLMM_High_Moderate.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

#glmm high vs control
HvC_GLMM <- read.table('GLMM_Results_HighvCont.txt', header =TRUE)
alpha8 <- quantile(read.table('Permutations_GLMM/Perms_GLMM_HighvCont.txt')[,1], .005)
Sig_HvC_GLMM <- subset(HvC_GLMM, pval <= alpha8)
write.table(Sig_HvC_GLMM,"Sig_SNPs_GLMM_High_Control.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

#glm control vs mod
MvC_GLMM <- read.table('GLMM_Results_ModvCont.txt', header =TRUE)
alpha9 <- quantile(read.table('Permutations_GLMM/Perms_GLMM_ModvCont.txt')[,1], .005)
Sig_MvC_GLMM <- subset(MvC_GLMM, pval <= alpha9)
write.table(Sig_MvC_GLMM,"Sig_SNPs_GLMM_Moderate_Control.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

