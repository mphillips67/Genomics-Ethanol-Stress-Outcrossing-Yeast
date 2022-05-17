setwd("/nfs3/IB/Burke_Lab/GLMM_perms_ethanol/Mod_Cont_perms")
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

###moderate vs control shuffle###
maf <- paste("alt_ETH_rep",c(control,moderate),"_wk15",sep="")
cov <- paste("N_ETH_rep",c(control,moderate),"_wk15",sep="")
pops <- c(rbind(maf,cov))
moderate_control <- data[,c("chr","pos",pops)]
sum_maf <- rowSums(moderate_control[,seq(3,ncol(moderate_control),by=2)])
moderate_control <- cbind(sum_maf, moderate_control)
moderate_control <- subset(moderate_control, sum_maf > 0 & sum_maf < 40)
moderate_control <- moderate_control[,2:ncol(moderate_control)]
rm(sum_maf)
rm(data)

shuf_maf <- sample(names( moderate_control[,seq(3,ncol(moderate_control),by=2)]))
shuf_cov <- gsub("alt", "N", shuf_maf)
shuf <- c(rbind(shuf_maf,shuf_cov))
Shuf_MC <- moderate_control[,c("chr","pos",shuf)]
rm(moderate_control)

### prep for glmm
A0 <- Shuf_MC[,seq(3,ncol(Shuf_MC),by=2)] * Shuf_MC[,seq(4,ncol(Shuf_MC),by=2)]
a0 <- Shuf_MC[,seq(4,ncol(Shuf_MC),by=2)] - A0
pop <- names(A0)
treatment <- c(rep("control",20), rep("moderate",20))

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

out.file <- paste("GLMM_perm_MvC", args[1], ".txt", sep="")
write.table(smallest_p, out.file, quote= FALSE, row.names= FALSE,col.names= FALSE, sep="\t")




