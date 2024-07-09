setwd("~/Dropbox/Ethanol expression project/Genomic_analysis/Haplotypes/")

library(tidyverse)
library(limSolve)
blah = function(par,tw){
  ww = exp(-unlist(tw)/(2*par^2))
  ww = ww/sum(ww)
  # the idea is to choose sigma, such that the 30 closest sites to pos
  # account for 75% of the weight
  (sum(sort(ww,decreasing=TRUE)[50:length(ww)])-0.50)^2
}



runscan <- function(inputfile, chr, poolroot, founderfile, outdir){
  # inputfile = file of SNP frequencies in pool and founders.  The structure of the file matters
  #     the 1st column is a chromosome name, and the 2nd a position, there can be additional columns
  #     eventually you want paired columns of the 1st is the frequency in some sample and the 2nd the coverage
  #     I have a python script that makes these files from "vcf" files, vcf are a very standard format
  #     this file should be preprocessed with GrannySNP.R
  # chr = chromosome to estimate frequencies for, e.g., "chrXI"
  # poolnumber = column number of frequency estimates for pool of interest
  # a file with founder information.  number of rows = number of founders.  First column is name in bam file (= root in SNP table), 2nd column is a short name
  # outdir = a prefix for the output file to write results to
  #     this greatly simplifies your life if you run this script for differents samples and chromsomes on a cluster...
  
  # inputfile = file of SNP frequencies in pool and founders.  The structure of the file matters
  #     the 1st column is a chromosome name, and the 2nd a position, there can be additional columns
  #     eventually you want paired columns of the 1st is the frequency in some sample and the 2nd the coverage
  #     I have a python script that makes these files from "vcf" files, vcf are a very standard format
  #     this file should be preprocessed with GrannySNP.R
  # chr = chromosome to estimate frequencies for, e.g., "chrXI"
  # poolnumber = column number of frequency estimates for pool of interest
  # a file with founder information.  number of rows = number of founders.  First column is name in bam file (= root in SNP table), 2nd column is a short name
  # outdir = a prefix for the output file to write results to
  #     this greatly simplifies your life if you run this script for differents samples and chromsomes on a cluster...
  
  yy = read.table(file=founderfile,header=FALSE)
  founderroots = yy[,1]
  foundershortnames = yy[,2]
  filename=paste(outdir,"/",poolroot,"_",chr,"_hap_freq.txt",sep='')
  
  if (!file.exists(filename)){	
    
    xx = read.table(file=inputfile,header=TRUE)
    pool <-  match(poolroot,names(xx))
    founders <- match(founderroots,names(xx))
    Nf = length(founders)
    # drop SNPs that are polymorphic in founders
    seg = apply(xx[,founders],1,function(x) sum(x>0.05 & x<0.95))
    xx = xx[seg==0,]
    
    seg = apply(xx[,founders],1,function(x) sum(unlist(lapply(x,function(y) min((1-y)^2,(0-y)^2)))))
    xx = xx[!is.na(seg),]
    seg = apply(xx[,founders],1,function(x) sum(unlist(lapply(x,function(y) min((1-y)^2,(0-y)^2)))))
    # this throws out about 2%
    xx = xx[seg<0.001,]
    
    #  here is something to try.
    # apparently this did not help
    # temp = round(xx[,founders],0)
    # xx[,founders] <- temp
    
    
    stepSize = 1000
    WINSIZE = 30000
    names = paste(foundershortnames, collapse=";")
    Maxpos = max(xx$POS[xx$CHROM == chr],na.rm=TRUE)
    Minpos = min(xx$POS[xx$CHROM == chr],na.rm=TRUE)
    ppp = seq(round((Minpos+5000)/100,0)*100,round((Maxpos-5000)/100,0)*100,stepSize)
    
    #remove bad window
    #ppp <- ppp[ppp!=941800 & ppp!=942800 & ppp!= 1135800  & ppp!=  1136800]
    
    lppp = length(ppp)
    ddd = data.frame(poolroot=rep("",lppp),chr=rep("",lppp),pos=rep(0,lppp),NSNPs=rep(0,lppp),foundernames=rep("",lppp),cutree=rep("",lppp),numberlevs=rep(0,lppp),founderfreqs=rep("",lppp),adjfounderfreqs=rep("",lppp),stringsAsFactors=FALSE)
    i=1
    
    
    
    for(pos in ppp){
      # print(pos)  #this is to check where bombs out
      
      #pos=9700
      #pos=50700
      # median non-recombined block size is 14kb
      #  so step size should be much smaller than that (500bp = 1/30th)
      #  I would guess optimum window size is +/- 7kb (= 14 total) ish.  
      #  but markers are limiting...
      predictors <- subset(xx, CHROM == chr & POS > pos - WINSIZE & POS < pos + WINSIZE, select=founders)
      Y <- subset(xx, CHROM == chr & POS > pos - WINSIZE & POS < pos + WINSIZE, select=pool)
      tw = (subset(xx, CHROM == chr & POS > pos - WINSIZE & POS < pos + WINSIZE, select=POS)-pos)^2
      # pick a sigma, near the ends of chromosomes sigma has to be bigger than in the middle
      predictNotMissing = apply(predictors,1,function(x) sum(is.na(x))==0) & (!is.na(Y))
      predictors = predictors[predictNotMissing,]
      Y = Y[predictNotMissing,]
      tw = tw[predictNotMissing,]
      out <- optimize(blah,c(2000,30000),tw)
      sigma=out$minimum
      weights = exp(-tw/(2*sigma^2))
      
      # clearly this is very much a heuristic
      # the current code basically will not calculate haplotype frequencies if two founders are really
      # similar in some window.  Another strategy would be to calculate haplotype frequencies, but print 
      # minManhattan.  Even better would be to think about which haplotypes are not distinguishable, and print that
      
      # Note that the haplotype frequencies are constrained to sum to one and have a minimum 
      # frequency of 0.002.  This means the maximum haplotpe frequency is 96.8%.
      # There are 17 haplotypes, if 16 have a frequency of 0.2% then the 17th
      # is 96.8%
      
      #test correlation matrix for na
      test <- cor(predictors*unlist(weights))
      test2 <- is.na(test) %>% table()
      
      MynumberSNP<- nrow(predictors)
      if(nrow(predictors)<10){
        myadjrsquare <- NA
        Freqs <- rep(NA,Nf)
        GroupAvFreqs <- rep(NA,Nf)
        Groups <- rep(NA,Nf)
      }else if (length(test2)==1){
        Groups <- cutree(hclust(as.dist(1-cor(predictors*unlist(weights))^2)),h=0.005)
        d = ncol(predictors)
        A = predictors
        B = Y
        E = t(matrix(rep(1,d)))
        F = 1
        G = diag(rep(1,d))
        H = matrix(rep(0.0003,d))
        #  G = diag(rep(1,d))
        #  H = matrix(rep(0,d))
        Wa = weights
        out = lsei(A=A,B=B,E=E,F=F,G=G,H=H,Wa=Wa,verbose=TRUE)
        myadjrsquare <- out$residualNorm                 
        Freqs <- out$X
        meanbygroup <- tapply(Freqs,as.numeric(Groups),mean)
        GroupAvFreqs <- meanbygroup[as.numeric(Groups)]      # freqs replaced by group means
        
        
      } # nrow > 10 
      
      groups = paste(as.numeric(Groups),collapse=";")
      freqs = paste(round(as.numeric(Freqs),4),collapse=";")
      afreqs = paste(round(as.numeric(GroupAvFreqs),4),collapse=";")
      nGroups = nlevels(as.factor(Groups))
      ddd[i,] = c(poolroot,chr,pos,MynumberSNP,names,groups,nGroups,freqs,afreqs)
      i=i+1
    }  # over positions
    write.table(ddd,filename,row.names=FALSE,sep="\t",quote=FALSE)
  } 
} # function specific for chr and pool

#read in files and paramters

inputfilename <- "Input_hapcaller_ethanol.txt" #name of input file
founderfile <- "Founders_12.txt" #name of founder file

moderate <- c(25:30,32:33,35:38,41:48)
week <- c("_wk01","_wk07")
poolroot <- paste("alt_ETH_rep",rep(moderate,each=2),week,sep="")
outdir <- "~/Dropbox/Ethanol expression project/Genomic_analysis/Haplotypes/Mod_Eth_1_7/"
data <- read.table("Input_hapcaller_ethanol.txt",header=T)

chr <- unique(data$CHROM)

for (j in 1:length(poolroot)){
  for(i in 1:length(chr)){
    runscan(inputfilename, as.character(chr[i]), poolroot[j], founderfile, outdir)
  }
}

