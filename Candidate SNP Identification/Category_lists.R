library(data.table)
library(dplyr)


#read in different candidate lists
control_cand <-  data.table(read.table("Sig_SNPs_Control_1v15.txt", header = TRUE))[,1:2]
mod_cand <-  data.table(read.table("Sig_SNPs_Moderate_1v15.txt", header = TRUE))[,1:2]
high_cand <-  data.table(read.table("Sig_SNPs_High_1v15.txt", header = TRUE))[,1:2]


#high specific 
temp <- anti_join(high_cand, mod_cand) #remove moderate candidates from high
high_only <- anti_join(temp, control_cand) #remove control candidates from high
write.table(high_only,"High_only_sig_snps.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

#moderate only
temp <- anti_join(mod_cand, high_cand) #remove high candidates from mod
mod_only <- anti_join(temp, control_cand) #remove control candidates from mod
write.table(mod_only,"Mod_only_sig_snps.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

#control only
temp <- anti_join(control_cand, high_cand) #remove high candidates from control
control_only <- anti_join(temp, mod_cand) #remove moderate candidates from control 
write.table(control_only,"Control_only_sig_snps.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

#outcrossing
temp <- merge(high_cand, control_cand, by = c("chr","pos")) #shared between high and control
outcrossing <- merge(temp, mod_cand, by = c("chr","pos")) #shared between high, control, and moderate
write.table(outcrossing,"Sig_SNPs_Outcrossing.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

#general ethanol 
temp <- merge(high_cand, mod_cand, by = c("chr","pos")) #overlap between high and mod candidates
ethanol <- anti_join(temp,control_cand) #removes sites sig change in controls
write.table(ethanol,"Sig_SNPs_Shared_High_Mod_Overlap.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

#high and control only
temp <- anti_join(high_cand, mod_cand) #remove moderate candidates from high
high_control <- merge(temp, control_cand) #snps shared with control
write.table(high_control,"High_Control_shared_snps.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

#moderate and control only
temp <- anti_join(mod_cand, high_cand) #remove high candidates from mod
mod_control <- merge(temp, control_cand) #snps shared with control
write.table(mod_control,"Mod_Control_shared_snps.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")





