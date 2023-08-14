setwd("~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/Modifed_CMH")
library(ACER)
library(poolSeq)

#read in SNP table
data <- read.table("~/Dropbox/Ethanol expression project/Genomic_analysis/Prep_quality_checks/SNPtable_raw.txt", header = TRUE)

#time points and reps
control <-c(paste("0",1:9,sep=""),10,12:21)
moderate <- c(25:30,32:33,35:38,41:48)
high <- c(49:56,59:61,64:72)
weeks <- c("_wk01","_wk07", "_wk15")


#general parameters that are the same for all comparisons (e.g pool size, reps, generations)
rep <- c(1:20)
tp <- c(1,7,15) * 14
ps <- rep(1e+7, 60)

#C populations
control_freq <- as.matrix(data[, paste("alt_ETH_rep", rep(control,each=3),weeks,sep="")]) #frequencies 
control_cov <- as.matrix(data[,paste("N_ETH_rep", rep(control,each=3),weeks,sep="")]) #coverage associated with frequencies
n_C <- round(c(261.91, 474.48, 331.84, 147.86, 422.94, 340.96, 432.57, 474.96, 253.29, 271.29, 435.57, 401.33, 389.45, 444.30, 397.15, 299.82, 433.46, 119.85, 283.84, 448.03))
pval <- adapted.cmh.test(freq=control_freq, coverage=control_cov, Ne= n_C, gen=tp, repl=rep, poolSize=ps, IntGen = TRUE)
control_results <- cbind(data[,1:2], pval)
write.table(control_results,"CMH_Results_Controls_1v15.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

#M populations
moderate_freq <- as.matrix(data[, paste("alt_ETH_rep", rep(moderate,each=3),weeks,sep="")]) #frequencies 
moderate_cov <- as.matrix(data[,paste("N_ETH_rep", rep(moderate,each=3),weeks,sep="")]) #coverage associated with frequencies
n_M <- round(c(542.78, 531.15, 488.03, 427.97, 496.75, 517.31, 559.11, 536.98, 532.67, 224.18, 305.01, 324.29, 568.15, 420.74, 338.30, 416.75, 515.01, 335.10, 480.61, 530.43))
pval <- adapted.cmh.test(freq=moderate_freq, coverage=moderate_cov, Ne= n_M, gen=tp, repl=rep, poolSize=ps, IntGen = TRUE)
moderate_results <- cbind(data[,1:2], pval)
write.table(moderate_results,"CMH_Results_Moderate_1v15.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

#H populations
high_freq <- as.matrix(data[, paste("alt_ETH_rep", rep(high,each=3),weeks,sep="")]) #frequencies 
high_cov <- as.matrix(data[,paste("N_ETH_rep", rep(high,each=3),weeks,sep="")]) #coverage associated with frequencies
n_H <- round(c(461.85, 430.33, 479.10, 542.53, 400.88, 363.14, 351.13, 467.39, 479.26, 469.96, 580.44, 523.56, 489.42, 513.20, 492.49, 439.40, 536.02, 585.46, 388.62, 434.82))
pval <- adapted.cmh.test(freq=high_freq, coverage=high_cov, Ne= n_H, gen=tp, repl=rep, poolSize=ps, IntGen = TRUE)
high_results <- cbind(data[,1:2], pval)
write.table(high_results,"CMH_Results_High_1v15.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")



