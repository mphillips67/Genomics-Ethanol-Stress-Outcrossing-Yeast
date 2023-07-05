setwd("~/Dropbox/Ethanol expression project/Genomic_analysis/Ne_S/Bait-ER")
library(cowplot)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(grid)
library(ggcorrplot)
library(dplyr)
library(GGally)

#read in results and make S table
control_s <- read.table("Control_BaitER_results.txt", header = TRUE)
mod_s <- read.table("Moderate_BaitER_results.txt", header = TRUE)
high_s <- read.table("High_BaitER_results.txt", header = TRUE)
data <- data.frame(control_s[,1:4], mod_s$sigma, high_s$sigma)
names(data) <- c("CHR", "BP", "REF", "C-populations","M-populations", "H-populations")


#read in cmh results
high_ethanol_cmh <- read.table('~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/CMH/CMH_Results_High_1v15.txt', header = TRUE, col.names = c("CHR","BP","P"))
mod_ethanol_cmh <- read.table('~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/CMH/CMH_Results_Moderate_1v15.txt', header = TRUE, col.names = c("CHR","BP","P"))
control_cmh <- read.table('~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/CMH/CMH_Results_Controls_1v15.txt', header = TRUE, col.names = c("CHR","BP","P"))

#read in category lists
high_only <- read.table("~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/Ethanol SNPs/High_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
mod_only <- read.table("~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/Ethanol SNPs/Mod_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
control_only <- read.table("~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/Ethanol SNPs/Control_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
general_lab  <- read.table("~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/Outcrossing Snps/Sig_SNPs_Outcrossing.txt", header =TRUE,col.names = c("CHR","BP"))
general_eth <- read.table("~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/Ethanol SNPs/Sig_SNPs_Shared_High_Mod_Overlap.txt", header =TRUE,col.names = c("CHR","BP"))


#genome wide s r2 plot
corr_all <- cor(data[,4:6], use = "complete.obs", method = "pearson") #cor matrix
all <- data[,4:6]
a <- ggpairs(all, title = "A. All Polymorphic Sites")


#general lab 
general_lab <- merge(general_lab, control_cmh, by = c("CHR","BP"))
MSM_lab <- general_lab %>%
  group_by(CHR) %>%
  filter(P == min(P)) %>%
  select(CHR, BP) #find MSM for each chr
extract_positions <- function(position, range) {
  filter(general_lab, CHR == position$CHR & BP >= position$BP - range & BP <= position$BP + range)
} #function to get sits in window around msm
range <- 2500 #window range - before and after
general_lab_msm_range <-bind_rows(lapply(seq_len(nrow(MSM_lab)), function(i) extract_positions(MSM_lab[i, ], range = range)))[,1:2]
lab <- na.omit(merge(data, general_lab_msm_range))[,4:6]
b <- ggpairs(lab, title = "B. Top General Lab Selection Candidates")

#general ethanol
general_eth <- merge(general_eth, high_ethanol_cmh, by = c("CHR","BP"))
MSM_eth <- general_eth %>%
  group_by(CHR) %>%
  filter(P == min(P)) %>%
  select(CHR, BP) #find MSM for each chr
extract_positions <- function(position, range) {
  filter(general_eth, CHR == position$CHR & BP >= position$BP - range & BP <= position$BP + range)
} #function to get sits in window around msm
range <- 2500 #window range - before and after
general_eth_msm_range <-bind_rows(lapply(seq_len(nrow(MSM_eth)), function(i) extract_positions(MSM_eth[i, ], range = range)))[,1:2]
eth <- na.omit(merge(data, general_eth_msm_range))[,4:6]
c <- ggpairs(eth, title = "C. Top General Ethanol Selection Candidates")

#high ethanol specific 
high_only <- merge(high_only, high_ethanol_cmh, by = c("CHR","BP"))
MSM_high <- high_only %>%
  group_by(CHR) %>%
  filter(P == min(P)) %>%
  select(CHR, BP) #find MSM for each chr
extract_positions <- function(position, range) {
  filter(high_only, CHR == position$CHR & BP >= position$BP - range & BP <= position$BP + range)
} #function to get sits in window around msm
range <- 2500 #window range - before and after
high_eth_msm_range <-bind_rows(lapply(seq_len(nrow(MSM_high)), function(i) extract_positions(MSM_high[i, ], range = range)))[,1:2]
high <- na.omit(merge(data, high_eth_msm_range))[,4:6]
d <- ggpairs(high, title = "D. Top High Ethanol Specific Candidates")

#moderate specific
mod_only <- merge(mod_only, mod_ethanol_cmh, by = c("CHR","BP"))
MSM_mod <- mod_only %>%
  group_by(CHR) %>%
  filter(P == min(P)) %>%
  select(CHR, BP) #find MSM for each chr
extract_positions <- function(position, range) {
  filter(mod_only, CHR == position$CHR & BP >= position$BP - range & BP <= position$BP + range)
} #function to get sits in window around msm
range <- 2500 #window range - before and after
mod_eth_msm_range <-bind_rows(lapply(seq_len(nrow(MSM_mod)), function(i) extract_positions(MSM_mod[i, ], range = range)))[,1:2]
mod <- na.omit(merge(data, mod_eth_msm_range))[,4:6]
e <- ggpairs(mod, title = "E. Top Moderate Ethanol Specific Candidates")

#control
control_only <- merge(control_only, control_cmh, by = c("CHR","BP"))
MSM_control <- control_only %>%
  group_by(CHR) %>%
  filter(P == min(P)) %>%
  select(CHR, BP) #find MSM for each chr
extract_positions <- function(position, range) {
  filter(control_only, CHR == position$CHR & BP >= position$BP - range & BP <= position$BP + range)
} #function to get sits in window around msm
range <- 2500 #window range - before and after
control_eth_msm_range <-bind_rows(lapply(seq_len(nrow(MSM_control)), function(i) extract_positions(MSM_control[i, ], range = range)))[,1:2]
control <- na.omit(merge(data, control_eth_msm_range))[,4:6]

f <- ggpairs(control, title = "F. Top Control Specific Candidates")


plot_grid(
  ggmatrix_gtable(a),
  ggmatrix_gtable(b),
  ggmatrix_gtable(c),
  ggmatrix_gtable(d),
  ggmatrix_gtable(e),
  ggmatrix_gtable(f),
  nrow = 2)
