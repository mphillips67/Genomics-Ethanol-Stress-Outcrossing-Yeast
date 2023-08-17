#script to generate average Ne estimates for each treatment. 
library(poolSeq)

data <- read.table("SNPtable_raw.txt", header = TRUE)

#ancestor data
anc_freq <- data$alt_ETH_anc_12SH_12
anc_cov <- data$N_ETH_anc_12SH_12
#treatment + tps
control <-c(paste("0",1:9,sep=""),10,12:21)
moderate <- c(25:30,32:33,35:38,41:48)
high <- c(49:56,59:61,64:72)
weeks <- c( "_wk15")

#control Ne
maf <- paste("alt_ETH_rep", rep(control,each=1),weeks,sep="")
cov <- paste("N_ETH_rep", rep(control,each=1),weeks,sep="")
pops <- c(rbind(maf,cov))
controls <- data[,c(pops)]
control_ne <- NULL
for (i in seq(1,ncol(controls),by=2)){
  
  z <- i+1
  temp <- controls[,i:z]
  ne <- estimateNe(p0 = anc_freq, pt = temp[,1], cov0 =  anc_cov, covt = temp[,2], t = 200,  ploidy = 2, method = "P.planII", poolSize = c(1e+7, 1e+7))
  control_ne <- rbind(control_ne,ne)
}
mean(control_ne)

#moderate Ne
maf <- paste("alt_ETH_rep", rep(moderate,each=1),weeks,sep="")
cov <- paste("N_ETH_rep", rep(moderate,each=1),weeks,sep="")
pops <- c(rbind(maf,cov))
mods<- data[,c(pops)]
mod_ne <- NULL
for (i in seq(1,ncol(mods),by=2)){
  
  z <- i+1
  temp <- mods[,i:z]
  ne <- estimateNe(p0 = anc_freq, pt = temp[,1], cov0 =  anc_cov, covt = temp[,2], t = 200,  ploidy = 2, method = "P.planII", poolSize = c(1e+7, 1e+7))
  mod_ne <- rbind(mod_ne,ne)
}
mean(mod_ne)

#high Ne
maf <- paste("alt_ETH_rep", rep(high,each=1),weeks,sep="")
cov <- paste("N_ETH_rep", rep(high,each=1),weeks,sep="")
pops <- c(rbind(maf,cov))
highs <- data[,c(pops)]
high_ne <- NULL
for (i in seq(1,ncol(highs),by=2)){
  
  z <- i+1
  temp <- highs[,i:z]
  ne <- estimateNe(p0 = anc_freq, pt = temp[,1], cov0 =  anc_cov, covt = temp[,2], t = 200,  ploidy = 2, method = "P.planII", poolSize = c(1e+7, 1e+7))
  high_ne <- rbind(high_ne,ne)
}
mean(high_ne)
