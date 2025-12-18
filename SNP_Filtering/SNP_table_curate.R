setwd("~/Dropbox/Ethanol expression project/Genomic_analysis/Prep_quality_checks/")
library("plyr")

data <- read.table("~/Dropbox/Ethanol expression project/Ethanol_snptables_091921/filtered_snps.txt", header = TRUE)


#filter based on min coverage in experimental populations
snpmin<-apply(data[,c(seq(53,ncol(data),by=2))],1,min)
temp <- cbind(data,snpmin)
SNP_cov_filt <- subset(temp, snpmin >= 20)

#filter based on possible SNPs given founders 
data2 <-  subset(data, CHROM != "chrmito")
x <- c("A1","A2","A6","A7","A8","A12","B3", "B4","B5","B7","B10","B11")
alt <- paste("alt_ETH_hap_", x, "_00", sep="")
cov <- paste("N_ETH_hap_", x, "_00", sep="")
haps <-c(rbind(alt,cov))
data3 <- data2[, c("CHROM", "POS", "REF", "ALT", haps)]
####calc maf
data3[,c(seq(5,ncol(data3)-1,by=2))] <- data3[,c(seq(5,ncol(data3)-1,by=2))]/data3[,c(seq(6,ncol(data3),by=2))]

#####filter things not fixed
for (i in seq(5,27,by=2)){
  data3 <- subset(data3, data3[i] == 0 |data3[i] == 1 )
} 
####make sure polymorphic (sum maf not 0 or 12)
polym <-apply(data3[,c(seq(5,ncol(data3),by=2))],1,sum)
data3 <- cbind(data3,polym)
data3 <- subset(data3, polym > 0)
data3 <- subset(data3, polym < 12)
possible <- data3[,1:2]
possible$id <- 1:nrow(possible)
SNP_cov_filt2 <- merge(SNP_cov_filt,possible, by = c("CHROM","POS"))
SNP_cov_filt2 <- SNP_cov_filt2[order(SNP_cov_filt2$id),]

####remove founder data and id column
#####treatment + tps
control <-c(paste("0",1:9,sep=""),10,12:21)
moderate <- c(25:30,32:33,35:38,41:48)
high <- c(49:56,59:61,64:72)
all_treatments <- c(control,moderate,high)
weeks <- c("_wk01","_wk07","_wk15")
maf <- paste("alt_ETH_rep", rep(all_treatments,each=3),weeks,sep="")
cov <- paste("N_ETH_rep", rep(all_treatments,each=3),weeks,sep="")
pops <- c(rbind(maf,cov))

SNP_cov_filt3 <- SNP_cov_filt2[,c("Nmiss","CHROM","POS","REF","ALT","alt_ETH_anc_12S_00","N_ETH_anc_12S_00","alt_ETH_anc_12S_06","N_ETH_anc_12S_06","alt_ETH_anc_12S_12", "N_ETH_anc_12S_12", "alt_ETH_anc_12SH_12","N_ETH_anc_12SH_12",pops,"snpmin")]


#filter cases where alt nuc is fixed and site not actually poly 
minor <- rowSums(SNP_cov_filt3[,c(seq(6,ncol(SNP_cov_filt3)-1,by=2))])/rowSums(SNP_cov_filt3[,c(seq(7,ncol(SNP_cov_filt3),by=2))]) #sum of alt/sum total cov
temp2 <- cbind(SNP_cov_filt3, minor)
Table_no_fixed_alt <- subset(temp2, minor != 1)

#Minor allele frequencey filter (2% "rule").
SNP_table_filt <- subset(Table_no_fixed_alt, minor >= 0.02 & minor <= 0.98)
SNP_table_filt <- SNP_table_filt[,2:373]

#alt nuc freq table
SNP_table_filt[,c(seq(5,ncol(SNP_table_filt)-1,by=2))] <- SNP_table_filt[,c(seq(5,ncol(SNP_table_filt)-1,by=2))]/SNP_table_filt[,c(seq(6,ncol(SNP_table_filt),by=2))]

#ancestor not fixed
SNP_table_filt <- subset(SNP_table_filt, alt_ETH_anc_12SH_12 > 0 & alt_ETH_anc_12SH_12 < 1)

#renaming chromosomes
SNP_table_filt$CHROM <- as.character(SNP_table_filt$CHROM)
replacement <- c('C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16')
old <- as.vector(unique(SNP_table_filt$CHROM))
for (i in 1:16){
  SNP_table_filt[SNP_table_filt==old[i]] <- replacement[i]
}


##rename first 4 cols 
names(SNP_table_filt)[1:4] <- c("chr", "pos", "ref", "alt")

#remove rows with unsual count pattern where present in ancestor, but missing in some cycle 1 reps despite being prsent at later gens
remove_chr <- "C04"
remove_pos <- c(540228, 520955, 520956)
SNP_table_filt <- subset(SNP_table_filt, !(chr == remove_chr & pos %in% remove_pos))


write.table(SNP_table_filt,"SNPtable_raw.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")









#scale to 50x coverage
temp <- SNP_table_filt
temp[,c(seq(6,ncol(temp),by=2))] <- 50
temp[,c(seq(5,ncol(temp),by=2))] <- round(temp[,c(seq(5,ncol(temp),by=2))]*temp[,c(seq(6,ncol(temp),by=2))])
temp[,c(seq(5,ncol(temp),by=2))] <- temp[,c(seq(5,ncol(temp),by=2))]/temp[,c(seq(6,ncol(temp),by=2))]

write.table(snp.data,"SNPtable_scales.txt", quote= FALSE, row.names= FALSE,col.names= TRUE, sep="\t")




#afs
#quick distribution plots
op <- par(mfrow = c(5,3),
          oma = c(5,4,0,0) + 0.1,
          mar = c(1.5,5,2,2) + 0.1)

alt <-  temp[,c(seq(5,ncol(temp),by=2))]

for (col in 1:ncol(alt)) {
  hist(alt[,col], main = names(alt[col]))
}


#genome het
het <- alt
for (i in 1:ncol(alt)){
  het[,i] <- 2 * het[,i] * (1-het[,i])
}


#pca
data <- t(as.matrix(alt))
pca <- prcomp(data, scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

Population <- c("Anc0",rep(c("Y1","Y7","Y15"),times=20),rep(c("M1","M7","M15"),times=20),rep(c("H1","H7","H15"),times=20))

library(ggplot2)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data
####labled by pop
ggplot(data=pca.data, aes(x=X, y=Y, label=Population, color = Population)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("")



