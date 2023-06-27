
#read in table + make table with only alt freqs
data <- read.table("SNPtable_50X_scaled.txt",header= TRUE)
alt <-  data[,c(seq(11,ncol(data),by=2))]

#pca plot
data <- t(as.matrix(alt))
pca <- prcomp(data, scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

Population <- c("Ancestor",rep(c("Control Ethanol C1","Control Ethanol C7","Control Ethanol C15"),times=20),rep(c("Moderate Ethanol C1","Moderate Ethanol C7","Moderate Ethanol C15"),times=20),rep(c("High Ethanol C1","High Ethanol C7","High Ethanol C15"),times=20))

library(ggplot2)
library(ggh4x)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])


p<-ggplot(pca.data,aes(x=X, y=Y,color=Population, label=Population ))
  p<-p+geom_point() + xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) + ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) + theme(legend.position = "none") + stat_centroid(aes(label = Population),size=6,check_overlap=T,geom = "text") 
p
#+ scale_color_manual(values=c("#56B4E9","#009E73","#D55E00")) #+stat_ellipse()
p