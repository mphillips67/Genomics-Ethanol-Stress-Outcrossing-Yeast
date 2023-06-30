#requires cmh, s, and category files. makes plots showing s estimate correlations between groups of experimental populations for top candiates in each response category. top candidates defined as sites within 5KB window around most singificant marker for largeest peak on each chromsome across the genome. 

library(gridExtra)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(grid)
library(ggcorrplot)
library(dplyr)

#read in results and make S table
control_s <- read.table("Control_BaitER_results.txt", header = TRUE)
mod_s <- read.table("Moderate_BaitER_results.txt", header = TRUE)
high_s <- read.table("High_BaitER_results.txt", header = TRUE)
data <- data.frame(control_s[,1:4], mod_s$sigma, high_s$sigma)
names(data) <- c("CHR", "BP", "REF", "C-populations","M-populations", "H-populations")


#read in cmh results
high_ethanol_cmh <- read.table('CMH_Results_High_1v15.txt', header = TRUE, col.names = c("CHR","BP","P"))
mod_ethanol_cmh <- read.table('CMH_Results_Moderate_1v15.txt', header = TRUE, col.names = c("CHR","BP","P"))
control_cmh <- read.table('CMH_Results_Controls_1v15.txt', header = TRUE, col.names = c("CHR","BP","P"))

#read in category lists
high_only <- read.table("High_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
mod_only <- read.table("Mod_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
control_only <- read.table("Control_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
general_lab  <- read.table("Sig_SNPs_Outcrossing.txt", header =TRUE,col.names = c("CHR","BP"))
general_eth <- read.table("Sig_SNPs_Shared_High_Mod_Overlap.txt", header =TRUE,col.names = c("CHR","BP"))


#genome wide s r2 plot
corr_all <- cor(data[,4:6], use = "complete.obs", method = "pearson") #cor matrix
all <- data[,4:6]
a <- ggcorrplot(corr_all, digits= 2, method = "square", colors = c("#E46726", "white", "#6D9EC1"), outline.col = "white", lab = TRUE, type = "lower", title = "A. All Polymorphic Sites") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))


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
corr_lab <- cor(lab, use = "complete.obs", method = "pearson")

b <- ggcorrplot(corr_lab, digits= 2,method = "square", colors = c("#E46726", "white", "#6D9EC1"), outline.col = "white", lab = TRUE, type = "lower",title = "B. Top General Lab Selection Candidates") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

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
corr_eth <- cor(eth, use = "complete.obs", method = "pearson")

c <- ggcorrplot(corr_eth, digits= 2,method = "square", colors = c("#E46726", "white", "#6D9EC1"), outline.col = "white", lab = TRUE, type = "lower",title = "C. Top General Ethanol Selection Candidates") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

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
corr_high <- cor(high, use = "complete.obs", method  = "pearson")

d <- ggcorrplot(corr_high, digits= 2,method = "square", colors = c("#E46726", "white", "#6D9EC1"), outline.col = "white", lab = TRUE, type = "lower",title = "D. Top High Ethanol Specific Candidates") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

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
corr_mod <- cor(mod, use = "complete.obs", method  = "pearson")

e <- ggcorrplot(corr_mod, digits= 2,method = "square", colors = c("#E46726", "white", "#6D9EC1"), outline.col = "white", lab = TRUE, type = "lower",title = "E. Top Moderate Ethanol Specific Candidates") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))


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
corr_control <- cor(control, use = "complete.obs", method  = "pearson")

f <- ggcorrplot(corr_control, digits= 2,method = "square", colors = c("#E46726", "white", "#6D9EC1"), outline.col = "white", lab = TRUE, type = "lower",title = "F. Top Control Specific Candidates") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

grid.arrange(a,b,c,d,e,f,nrow=2)
