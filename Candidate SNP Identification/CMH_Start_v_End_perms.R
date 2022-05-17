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
####cmh perms
pval_control <- NULL
for (i in 1:5000){
shuf_maf <- sample(names( control_cmh[,seq(3,ncol(control_cmh),by=2)]))
shuf_cov <- gsub("alt", "N", shuf_maf)
shuf <- c(rbind(shuf_maf,shuf_cov))
temp <- data[,c("chr","pos",shuf)]
A0 <- temp[,seq(3,ncol(temp),by=4)] * temp[,seq(4,ncol(temp),by=4)]
a0 <- temp[,seq(4,ncol(temp),by=4)] - A0
A0 <- t(A0)
a0 <- t(a0)
At <- temp[,seq(5,ncol(temp),by=4)] * temp[,seq(6,ncol(temp),by=4)]
at <- temp[,seq(6,ncol(temp),by=4)] - At
At <- t(At)
at <- t(at)
pval <- cmh.test(A0=A0, a0=a0, At=At, at=at)
temp2 <- na.omit(cbind(temp[,1:2], pval))
temp2 <- subset(temp2, chr != "chrmito")
pval_control <- rbind(pval_control,min(temp2$pval,na.rm = TRUE))
rm(temp,temp2,A0,a0,At,at)
}
write.table(pval_control,"Perms_Control_1v15_5k.txt", quote= FALSE, row.names= FALSE,col.names= FALSE, sep="\t")
rm(control_cmh)


#moderate start vs. end####
maf <- paste("alt_ETH_rep", rep(moderate,each=2),weeks,sep="")
cov <- paste("N_ETH_rep", rep(moderate,each=2),weeks,sep="")
pops <- c(rbind(maf,cov))
moderate_cmh <- data[,c("chr","pos",pops)]
####cmh perms
pval_moderate <- NULL
for (i in 1:5000){
  shuf_maf <- sample(names( moderate_cmh[,seq(3,ncol(moderate_cmh),by=2)]))
  shuf_cov <- gsub("alt", "N", shuf_maf)
  shuf <- c(rbind(shuf_maf,shuf_cov))
  temp <- data[,c("chr","pos",shuf)]
  A0 <- temp[,seq(3,ncol(temp),by=4)] * temp[,seq(4,ncol(temp),by=4)]
  a0 <- temp[,seq(4,ncol(temp),by=4)] - A0
  A0 <- t(A0)
  a0 <- t(a0)
  At <- temp[,seq(5,ncol(temp),by=4)] * temp[,seq(6,ncol(temp),by=4)]
  at <- temp[,seq(6,ncol(temp),by=4)] - At
  At <- t(At)
  at <- t(at)
  pval <- cmh.test(A0=A0, a0=a0, At=At, at=at)
  temp2 <- na.omit(cbind(temp[,1:2], pval))
  temp2 <- subset(temp2, chr != "chrmito")
  pval_moderate <- rbind(pval_moderate,min(temp2$pval,na.rm = TRUE))
  rm(temp,temp2,A0,a0,At,at)
}
write.table(pval_moderate,"Perms_Moderate_1v15_5k.txt", quote= FALSE, row.names= FALSE,col.names= FALSE, sep="\t")
rm(moderate_cmh)

###high start vs. end###
maf <- paste("alt_ETH_rep", rep(high,each=2),weeks,sep="")
cov <- paste("N_ETH_rep", rep(high,each=2),weeks,sep="")
pops <- c(rbind(maf,cov))
high_cmh <- data[,c("chr","pos",pops)]
pval_high <- NULL
for (i in 1:5000){
  shuf_maf <- sample(names( high_cmh[,seq(3,ncol(high_cmh),by=2)]))
  shuf_cov <- gsub("alt", "N", shuf_maf)
  shuf <- c(rbind(shuf_maf,shuf_cov))
  temp <- data[,c("chr","pos",shuf)]
  A0 <- temp[,seq(3,ncol(temp),by=4)] * temp[,seq(4,ncol(temp),by=4)]
  a0 <- temp[,seq(4,ncol(temp),by=4)] - A0
  A0 <- t(A0)
  a0 <- t(a0)
  At <- temp[,seq(5,ncol(temp),by=4)] * temp[,seq(6,ncol(temp),by=4)]
  at <- temp[,seq(6,ncol(temp),by=4)] - At
  At <- t(At)
  at <- t(at)
  pval <- cmh.test(A0=A0, a0=a0, At=At, at=at)
  temp2 <- na.omit(cbind(temp[,1:2], pval))
  temp2 <- subset(temp2, chr != "chrmito")
  pval_high <- rbind(pval_high,min(temp2$pval,na.rm = TRUE))
  rm(temp,temp2,A0,a0,At,at)
}
write.table(pval_high,"Perms_High_1v15_5k.txt", quote= FALSE, row.names= FALSE,col.names= FALSE, sep="\t")
rm(high_cmh)


