setwd("/nfs3/IB/Burke_Lab/GLMM_perms_ethanol/High_Mod_perms")
library(lme4)
library(plyr)
args <- commandArgs(TRUE)
print(args[1])

#read in table
data <- read.table("SNPtable_50X_scaled.txt", header = TRUE)

#treatments
control <-c(paste("0",1:9,sep=""),10,12:21)
moderate <- c(25:30,32:33,35:38,41:48)
high <- c(49:56,59:61,64:72)

###high vs moderate shuffle###
maf <- paste("alt_ETH_rep",c(moderate,high),"_wk15",sep="")
cov <- paste("N_ETH_rep",c(moderate,high),"_wk15",sep="")
pops <- c(rbind(maf,cov))
high_moderate <- data[,c("chr","pos",pops)]
sum_maf <- rowSums(high_moderate[,seq(3,ncol(high_moderate),by=2)])
high_moderate <- cbind(sum_maf, high_moderate)
high_moderate <- subset(high_moderate, sum_maf > 0 & sum_maf < 40)
high_moderate <- high_moderate[,2:ncol(high_moderate)]
rm(sum_maf)
rm(data)

shuf_maf <- sample(names( high_moderate[,seq(3,ncol(high_moderate),by=2)]))
shuf_cov <- gsub("alt", "N", shuf_maf)
shuf <- c(rbind(shuf_maf,shuf_cov))
Shuf_HM <- high_moderate[,c("chr","pos",shuf)]
rm(high_moderate)

### prep for glmm
A0 <- Shuf_HM[,seq(3,ncol(Shuf_HM),by=2)] * Shuf_HM[,seq(4,ncol(Shuf_HM),by=2)]
a0 <- Shuf_HM[,seq(4,ncol(Shuf_HM),by=2)] - A0
pop <- names(A0)
treatment <- c(rep("moderate",20), rep("high",20))

## run glmm and record smalled p-value
pval <- numeric (nrow(A0))
for (i in 1:nrow(A0)){
  A1 <- t(A0[i,])
  A2 <- t(a0[i,])
  test <- suppressMessages(glmer(cbind(A1,A2)~treatment+(1|pop),family="binomial"))
  test2 <- suppressMessages(glmer(cbind(A1,A2)~(1|pop),family="binomial"))
  test3 <- anova(test,test2)
  pval[i] <- test3$`Pr(>Chisq)`[2]
  #print(i)
  rm(A1,A2)
}
smallest_p <- min(pval)

out.file <- paste("GLMM_perm_HvM", args[1], ".txt", sep="")
write.table(smallest_p, out.file, quote= FALSE, row.names= FALSE,col.names= FALSE, sep="\t")

