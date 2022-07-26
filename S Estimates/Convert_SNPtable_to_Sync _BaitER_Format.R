
data <- read.table("SNPtable_raw.txt", header = TRUE)

#identify high coverage sites
snpmax <- apply(data[,c(seq(6,ncol(data),by=2))],1,max)
temp <- cbind(data,snpmax)
data <- subset(temp, snpmax <= 800)[,1:372]

#ancestor data
anc_freq <- rep("alt_ETH_anc_12SH_12",20)
anc_cov <- rep("N_ETH_anc_12SH_12",20) 
ancestor <- c(rbind(anc_freq,anc_cov))

#treatment + tps
control <-c(paste("0",1:9,sep=""),10,12:21)
moderate <- c(25:30,32:33,35:38,41:48)
high <- c(49:56,59:61,64:72)
weeks <- c("_wk01","_wk07", "_wk15")

#control sync file
maf <- paste("alt_ETH_rep", control,rep(weeks,each=20),sep="")
cov <- paste("N_ETH_rep", control,rep(weeks,each=20),sep="")
pops <- c(rbind(maf,cov))
controls <- data[,c("chr","pos","ref","alt",ancestor,pops)]

#### replace freq with alt count and N with ref
controls[,c(seq(5,ncol(controls),by=2))] <- round(controls[,c(seq(5,ncol(controls),by=2))] * controls[,c(seq(6,ncol(controls),by=2))])
controls[,c(seq(6,ncol(controls),by=2))] <- controls[,c(seq(6,ncol(controls),by=2))] - controls[,c(seq(5,ncol(controls),by=2))] 

####create sync
control_sync <- cbind(controls[,1:3],data.frame(matrix(ncol = 80, nrow = nrow(controls))))
d <- 3
for (j in seq(5,ncol(controls),by = 2)){
  t <- j + 1
  temp <- rep(NA, nrow(controls))
  
  for (i in 1:nrow(controls)){
     if (as.character(controls[i,3]) =="A" && as.character(controls[i,4])=="T"){
        temp2 <-  paste(as.character(controls[i,t]),as.character(controls[i,j]),"0","0","0","0", sep=":")
        temp[i] <- temp2
     } else if (as.character(controls[i,3]) =="A" && as.character(controls[i,4])=="C"){
       temp2 <-  paste(as.character(controls[i,t]),"0",as.character(controls[i,j]),"0","0","0", sep=":")
       temp[i] <- temp2
     }else if (as.character(controls[i,3]) =="A" && as.character(controls[i,4])=="G"){
       temp2 <-  paste(as.character(controls[i,t]),"0","0",as.character(controls[i,j]),"0","0", sep=":")
       temp[i] <- temp2
     }else if (as.character(controls[i,3]) =="T" && as.character(controls[i,4])=="A"){
       temp2 <-  paste(as.character(controls[i,j]), as.character(controls[i,t]),"0","0","0","0", sep=":")
       temp[i] <- temp2
     }else if (as.character(controls[i,3]) =="C" && as.character(controls[i,4])=="A"){
       temp2 <- paste(as.character(controls[i,j]),"0",as.character(controls[i,t]),"0","0","0", sep=":")
       temp[i] <- temp2
     }else if (as.character(controls[i,3]) =="G" && as.character(controls[i,4])=="A"){
       temp2 <-  paste(as.character(controls[i,j]), "0","0",as.character(controls[i,t]),"0","0", sep=":")
       temp[i] <- temp2
     }else if (as.character(controls[i,3]) =="T" && as.character(controls[i,4])=="C"){
       temp2 <-  paste("0", as.character(controls[i,t]),as.character(controls[i,j]),"0","0","0", sep=":")
       temp[i] <- temp2
     }else if (as.character(controls[i,3]) =="T" && as.character(controls[i,4])=="G"){
       temp2 <-  paste("0", as.character(controls[i,t]),"0",as.character(controls[i,j]),"0","0", sep=":")
       temp[i] <- temp2
     }else if (as.character(controls[i,3]) =="C" && as.character(controls[i,4])=="T"){
       temp2 <-  paste("0", as.character(controls[i,j]),as.character(controls[i,t]),"0","0","0", sep=":")
       temp[i] <- temp2
     }else if (as.character(controls[i,3]) =="G" && as.character(controls[i,4])=="T"){
       temp2 <-  paste("0", as.character(controls[i,j]),"0",as.character(controls[i,t]),"0","0", sep=":")
       temp[i] <- temp2
     }else if (as.character(controls[i,3]) =="C" && as.character(controls[i,4])=="G"){
       temp2 <-  paste("0", "0",as.character(controls[i,t]),as.character(controls[i,j]),"0","0", sep=":")
       temp[i] <- temp2
     }else if (as.character(controls[i,3]) =="G" && as.character(controls[i,4])=="C"){
       temp2 <-  paste("0", "0",as.character(controls[i,j]),as.character(controls[i,t]),"0","0", sep=":")
       temp[i] <- temp2
     }
  } 
  d <- d +1
  #print(d)
  control_sync[,d] <- temp
  names(control_sync)[d] <- colnames(controls)[t]
}
write.table (control_sync,file="Controls_ypd_BaitER.sync",quote=FALSE,row.names = FALSE, col.names = TRUE, sep = "\t")



#moderate sync
maf <- paste("alt_ETH_rep", moderate,rep(weeks,each=20),sep="")
cov <- paste("N_ETH_rep", moderate,rep(weeks,each=20),sep="")
pops <- c(rbind(maf,cov))
mods <- data[,c("chr","pos","ref","alt",ancestor,pops)]

#### replace freq with alt count and N with ref
mods[,c(seq(5,ncol(mods),by=2))] <- round(mods[,c(seq(5,ncol(mods),by=2))] * mods[,c(seq(6,ncol(mods),by=2))])
mods[,c(seq(6,ncol(mods),by=2))] <- mods[,c(seq(6,ncol(mods),by=2))] - mods[,c(seq(5,ncol(mods),by=2))] 

####create sync
mod_sync <- cbind(mods[,1:3],data.frame(matrix(ncol = 80, nrow = nrow(mods))))
d <- 3
for (j in seq(5,ncol(mods),by = 2)){
  t <- j + 1
  temp <- rep(NA, nrow(mods))
  
  for (i in 1:nrow(mods)){
    if (as.character(mods[i,3]) =="A" && as.character(mods[i,4])=="T"){
      temp2 <-  paste(as.character(mods[i,t]),as.character(mods[i,j]),"0","0","0","0", sep=":")
      temp[i] <- temp2
    } else if (as.character(mods[i,3]) =="A" && as.character(mods[i,4])=="C"){
      temp2 <-  paste(as.character(mods[i,t]),"0",as.character(mods[i,j]),"0","0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(mods[i,3]) =="A" && as.character(mods[i,4])=="G"){
      temp2 <-  paste(as.character(mods[i,t]),"0","0",as.character(mods[i,j]),"0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(mods[i,3]) =="T" && as.character(mods[i,4])=="A"){
      temp2 <-  paste(as.character(mods[i,j]), as.character(mods[i,t]),"0","0","0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(mods[i,3]) =="C" && as.character(mods[i,4])=="A"){
      temp2 <- paste(as.character(mods[i,j]),"0",as.character(mods[i,t]),"0","0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(mods[i,3]) =="G" && as.character(mods[i,4])=="A"){
      temp2 <-  paste(as.character(mods[i,j]), "0","0",as.character(mods[i,t]),"0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(mods[i,3]) =="T" && as.character(mods[i,4])=="C"){
      temp2 <-  paste("0", as.character(mods[i,t]),as.character(mods[i,j]),"0","0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(mods[i,3]) =="T" && as.character(mods[i,4])=="G"){
      temp2 <-  paste("0", as.character(mods[i,t]),"0",as.character(mods[i,j]),"0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(mods[i,3]) =="C" && as.character(mods[i,4])=="T"){
      temp2 <-  paste("0", as.character(mods[i,j]),as.character(mods[i,t]),"0","0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(mods[i,3]) =="G" && as.character(mods[i,4])=="T"){
      temp2 <-  paste("0", as.character(mods[i,j]),"0",as.character(mods[i,t]),"0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(mods[i,3]) =="C" && as.character(mods[i,4])=="G"){
      temp2 <-  paste("0", "0",as.character(mods[i,t]),as.character(mods[i,j]),"0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(mods[i,3]) =="G" && as.character(mods[i,4])=="C"){
      temp2 <-  paste("0", "0",as.character(mods[i,j]),as.character(mods[i,t]),"0","0", sep=":")
      temp[i] <- temp2
    }
  } 
  d <- d +1
  #print(d)
  mod_sync[,d] <- temp
  names(mod_sync)[d] <- colnames(mods)[t]
}
write.table (mod_sync,file="Mod_eth_BaitER.sync",quote=FALSE,row.names = FALSE, col.names = TRUE, sep = "\t")



#high sync
maf <- paste("alt_ETH_rep", high,rep(weeks,each=20),sep="")
cov <- paste("N_ETH_rep", high,rep(weeks,each=20),sep="")
pops <- c(rbind(maf,cov))
highs <- data[,c("chr","pos","ref","alt",ancestor,pops)]

#### replace freq with alt count and N with ref
highs[,c(seq(5,ncol(highs),by=2))] <- round(highs[,c(seq(5,ncol(highs),by=2))] * highs[,c(seq(6,ncol(highs),by=2))])
highs[,c(seq(6,ncol(highs),by=2))] <- highs[,c(seq(6,ncol(highs),by=2))] - highs[,c(seq(5,ncol(highs),by=2))] 

####create sync
high_sync <- cbind(highs[,1:3],data.frame(matrix(ncol = 80, nrow = nrow(highs))))
d <- 3
for (j in seq(5,ncol(highs),by = 2)){
  t <- j + 1
  temp <- rep(NA, nrow(highs))
  
  for (i in 1:nrow(highs)){
    if (as.character(highs[i,3]) =="A" && as.character(highs[i,4])=="T"){
      temp2 <-  paste(as.character(highs[i,t]),as.character(highs[i,j]),"0","0","0","0", sep=":")
      temp[i] <- temp2
    } else if (as.character(highs[i,3]) =="A" && as.character(highs[i,4])=="C"){
      temp2 <-  paste(as.character(highs[i,t]),"0",as.character(highs[i,j]),"0","0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(highs[i,3]) =="A" && as.character(highs[i,4])=="G"){
      temp2 <-  paste(as.character(highs[i,t]),"0","0",as.character(highs[i,j]),"0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(highs[i,3]) =="T" && as.character(highs[i,4])=="A"){
      temp2 <-  paste(as.character(highs[i,j]), as.character(highs[i,t]),"0","0","0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(highs[i,3]) =="C" && as.character(highs[i,4])=="A"){
      temp2 <- paste(as.character(highs[i,j]),"0",as.character(highs[i,t]),"0","0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(highs[i,3]) =="G" && as.character(highs[i,4])=="A"){
      temp2 <-  paste(as.character(highs[i,j]), "0","0",as.character(highs[i,t]),"0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(highs[i,3]) =="T" && as.character(highs[i,4])=="C"){
      temp2 <-  paste("0", as.character(highs[i,t]),as.character(highs[i,j]),"0","0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(highs[i,3]) =="T" && as.character(highs[i,4])=="G"){
      temp2 <-  paste("0", as.character(highs[i,t]),"0",as.character(highs[i,j]),"0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(highs[i,3]) =="C" && as.character(highs[i,4])=="T"){
      temp2 <-  paste("0", as.character(highs[i,j]),as.character(highs[i,t]),"0","0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(highs[i,3]) =="G" && as.character(highs[i,4])=="T"){
      temp2 <-  paste("0", as.character(highs[i,j]),"0",as.character(highs[i,t]),"0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(highs[i,3]) =="C" && as.character(highs[i,4])=="G"){
      temp2 <-  paste("0", "0",as.character(highs[i,t]),as.character(highs[i,j]),"0","0", sep=":")
      temp[i] <- temp2
    }else if (as.character(highs[i,3]) =="G" && as.character(highs[i,4])=="C"){
      temp2 <-  paste("0", "0",as.character(highs[i,j]),as.character(highs[i,t]),"0","0", sep=":")
      temp[i] <- temp2
    }
  } 
  d <- d +1
  #print(d)
  high_sync[,d] <- temp
  names(high_sync)[d] <- colnames(highs)[t]
}
write.table (high_sync,file="High_eth_BaitER.sync",quote=FALSE,row.names = FALSE, col.names = TRUE, sep = "\t")
