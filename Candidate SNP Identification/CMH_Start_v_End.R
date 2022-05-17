
library(poolSeq)
#read in table
data <- read.table("SNPtable_50X_scaled.txt", header = TRUE)

#treatment + tps
control <-c(paste("0",1:9,sep=""),10,12:21)
moderate <- c(25:30,32:33,35:38,41:48)
high <- c(49:56,59:61,64:72)
weeks <- c("_wk01", "_wk15")

#control start vs. end####
maf <- paste("alt_ETH_rep", rep(control,each=2),weeks,sep="")
cov <- paste("N_ETH_rep", rep(control,each=2),weeks,sep="")
pops <- c(rbind(maf,cov))
control_cmh <- data[,c("chr","pos",pops)]

####prep inputs and run cmh
A0 <- control_cmh[,seq(3,ncol(control_cmh),by=4)] * control_cmh[,seq(4,ncol(control_cmh),by=4)]
a0 <- control_cmh[,seq(4,ncol(control_cmh),by=4)] - A0
A0 <- t(A0)
a0 <- t(a0)
At <- control_cmh[,seq(5,ncol(control_cmh),by=4)] * control_cmh[,seq(6,ncol(control_cmh),by=4)]
at <- control_cmh[,seq(6,ncol(control_cmh),by=4)] - At
At <- t(At)
at <- t(at)
pval <- cmh.test(A0=A0, a0=a0, At=At, at=at)
cmh_results_control <- na.omit(cbind(control_cmh[,1:2], pval))
cmh_results_control <- subset(cmh_results_control, chr != "chrmito")
write.table(cmh_results_control,"CMH_Results_Controls_1v15.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

#moderate start vs. end####
maf <- paste("alt_ETH_rep", rep(moderate,each=2),weeks,sep="")
cov <- paste("N_ETH_rep", rep(moderate,each=2),weeks,sep="")
pops <- c(rbind(maf,cov))
moderate_cmh <- data[,c("chr","pos",pops)]

####prep inputs and run cmh
A0 <- moderate_cmh[,seq(3,ncol(moderate_cmh),by=4)] * moderate_cmh[,seq(4,ncol(moderate_cmh),by=4)]
a0 <- moderate_cmh[,seq(4,ncol(moderate_cmh),by=4)] - A0
A0 <- t(A0)
a0 <- t(a0)
At <- moderate_cmh[,seq(5,ncol(moderate_cmh),by=4)] * moderate_cmh[,seq(6,ncol(moderate_cmh),by=4)]
at <- moderate_cmh[,seq(6,ncol(moderate_cmh),by=4)] - At
At <- t(At)
at <- t(at)
pval <- cmh.test(A0=A0, a0=a0, At=At, at=at)
cmh_results_moderate <- na.omit(cbind(moderate_cmh[,1:2], pval))
cmh_results_moderate <- subset(cmh_results_moderate, chr != "chrmito")
write.table(cmh_results_moderate,"CMH_Results_Moderate_1v15.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")

###high start vs. end###
maf <- paste("alt_ETH_rep", rep(high,each=2),weeks,sep="")
cov <- paste("N_ETH_rep", rep(high,each=2),weeks,sep="")
pops <- c(rbind(maf,cov))
high_cmh <- data[,c("chr","pos",pops)]

####prep inputs and run cmh
A0 <- high_cmh[,seq(3,ncol(high_cmh),by=4)] * high_cmh[,seq(4,ncol(high_cmh),by=4)]
a0 <- high_cmh[,seq(4,ncol(high_cmh),by=4)] - A0
A0 <- t(A0)
a0 <- t(a0)
At <- high_cmh[,seq(5,ncol(high_cmh),by=4)] * high_cmh[,seq(6,ncol(high_cmh),by=4)]
at <- high_cmh[,seq(6,ncol(high_cmh),by=4)] - At
At <- t(At)
at <- t(at)
pval <- cmh.test(A0=A0, a0=a0, At=At, at=at)
cmh_results_high <- na.omit(cbind(high_cmh[,1:2], pval))
cmh_results_high <- subset(cmh_results_high, chr != "chrmito")
write.table(cmh_results_high,"CMH_Results_High_1v15.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")


