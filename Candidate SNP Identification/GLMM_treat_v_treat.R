
library(lme4)
library(plyr)
#read in table
data <- read.table("SNPtable_50X_scaled.txt", header = TRUE)

#add 1 to every 0
#for (k in seq(7,ncol(data),by=2)){
#for (j in 1:nrow(data)){
#  x <- data[j,k] 
#if (x==0){
#  data[j,k] <- data[j,k]+(1/50)
#}}}
#treatments
control <-c(paste("0",1:9,sep=""),10,12:21)
moderate <- c(25:30,32:33,35:38,41:48)
high <- c(49:56,59:61,64:72)

###high vs moderate###
maf <- paste("alt_ETH_rep",c(moderate,high),"_wk15",sep="")
cov <- paste("N_ETH_rep",c(moderate,high),"_wk15",sep="")
pops <- c(rbind(maf,cov))
high_moderate <- data[,c("chr","pos",pops)]
sum_maf <- rowSums(high_moderate[,seq(3,ncol(high_moderate),by=2)])
high_moderate <- cbind(sum_maf, high_moderate)
high_moderate <- subset(high_moderate, sum_maf > 0 & sum_maf < 40)
high_moderate <- high_moderate[,2:ncol(high_moderate)]

A0 <- high_moderate[,seq(3,ncol(high_moderate),by=2)] * high_moderate[,seq(4,ncol(high_moderate),by=2)]
a0 <- high_moderate[,seq(4,ncol(high_moderate),by=2)] - A0
pop <- names(A0)
treatment <- c(rep("moderate",20), rep("high",20))

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
high_mod_glm <- cbind(high_moderate[,1:2],pval)
write.table(high_mod_glm,"GLMM_Results_HighvMod.txt", quote= FALSE, row.names= FALSE,col.names=TRUE, sep="\t")


###high vs control###
maf <- paste("alt_ETH_rep",c(control,high),"_wk15",sep="")
cov <- paste("N_ETH_rep",c(control,high),"_wk15",sep="")
pops <- c(rbind(maf,cov))
high_control <- data[,c("chr","pos",pops)]
sum_maf <- rowSums(high_control[,seq(3,ncol(high_control),by=2)])
high_control <- cbind(sum_maf, high_control)
high_control <- subset(high_control, sum_maf > 0 & sum_maf < 40)
high_control <- high_control[,2:ncol(high_control)]

A0 <- high_control[,seq(3,ncol(high_control),by=2)] * high_control[,seq(4,ncol(high_control),by=2)]
a0 <- high_control[,seq(4,ncol(high_control),by=2)] - A0
pop <- names(A0)
treatment <- c(rep("control",20), rep("high",20))

pval <- numeric (nrow(A0))
for (i in 1:nrow(A0)){
  A1 <- t(A0[i,])
  A2 <- t(a0[i,])
  test <- suppressMessages(glmer(cbind(A1,A2)~treatment+(1|pop),family="binomial"))
  test2 <- suppressMessages(glmer(cbind(A1,A2)~(1|pop),family="binomial"))
  test3 <- anova(test,test2)
  pval[i] <- test3$`Pr(>Chisq)`[2]
  rm(A1,A2)
  #print(i)
}
high_cont_glm <- cbind(high_control[,1:2],pval)
write.table(high_cont_glm,"GLMM_Results_HighvCont.txt", quote= FALSE, row.names= FALSE,col.names=TRUE, sep="\t")

###moderate vs control###
maf <- paste("alt_ETH_rep",c(control,moderate),"_wk15",sep="")
cov <- paste("N_ETH_rep",c(control,moderate),"_wk15",sep="")
pops <- c(rbind(maf,cov))
moderate_control <- data[,c("chr","pos",pops)]
sum_maf <- rowSums(moderate_control[,seq(3,ncol(moderate_control),by=2)])
moderate_control <- cbind(sum_maf, moderate_control)
moderate_control <- subset(moderate_control, sum_maf > 0 & sum_maf < 40)
moderate_control <- moderate_control[,2:ncol(moderate_control)]

A0 <- moderate_control[,seq(3,ncol(moderate_control),by=2)] * moderate_control[,seq(4,ncol(moderate_control),by=2)]
a0 <- moderate_control[,seq(4,ncol(moderate_control),by=2)] - A0
pop <- names(A0)
treatment <- c(rep("control",20), rep("moderate",20))

pval <- numeric (nrow(A0))
for (i in 1:nrow(A0)){
  A1 <- t(A0[i,])
  A2 <- t(a0[i,])
  test <- suppressMessages(glmer(cbind(A1,A2)~treatment+(1|pop),family="binomial"))
  test2 <- suppressMessages(glmer(cbind(A1,A2)~(1|pop),family="binomial"))
  test3 <- anova(test,test2)
  pval[i] <- test3$`Pr(>Chisq)`[2]
  rm(A1,A2)
  #print(i)
}
mod_cont_glm <- cbind(moderate_control[,1:2],pval)
write.table(mod_cont_glm,"GLMM_Results_ModvCont.txt", quote= FALSE, row.names= FALSE,col.names=TRUE, sep="\t")





