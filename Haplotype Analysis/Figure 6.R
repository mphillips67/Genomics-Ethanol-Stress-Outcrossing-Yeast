setwd("~/Dropbox/Ethanol expression project/Genomic_analysis/Haplotypes/")
library(ggplot2)
library(ggpubr)
library(gginference)
library(dplyr)
library(gridExtra)
library(reshape2)
library(stringr)
library(tidyverse)
library(grid)

#candidate SNPs
control_sig <- read.table("~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/Modifed_CMH/Sig_SNPs_Control_1v15.txt", header = TRUE)
moderate_sig <- read.table("~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/Modifed_CMH/Sig_SNPs_Moderate_1v15.txt", header = TRUE)
high_sig <- read.table("~/Dropbox/Ethanol expression project/Genomic_analysis/SNP freq comparisons/Modifed_CMH/Sig_SNPs_High_1v15.txt" ,header = TRUE)

#general ethanol selection
#temp <- merge(high_sig, control_sig, by = c("chr","pos")) #shared between high and control
general <- merge(high_sig, moderate_sig, by = c("chr","pos")) #shared between high, control, and moderate
general <- subset(general, chr == "C02")
general_max <- general[which.min(general$pval.x), ]
general_max_pos <- general_max[,2]

#all haps 
control_c1 <-  read.table("Control//Control_wk1_haps.txt", header = TRUE)
founders <- as.vector(str_split_fixed(control_c1[1,5], ";",12))
control_c1_freqs <- cbind(control_c1[,1:3],colsplit(control_c1$adjfounderfreqs, ";", founders))
for (i in 4:ncol(control_c1_freqs)){
  control_c1_freqs <- subset(control_c1_freqs, control_c1_freqs[i] >= 0 & control_c1_freqs[i] <= 1 )
}
control_c7 <-  read.table("Control//Control_wk7_haps.txt", header = TRUE)
founders <- as.vector(str_split_fixed(control_c7[1,5], ";",12))
control_c7_freqs <- cbind(control_c7[,1:3],colsplit(control_c7$adjfounderfreqs, ";", founders))
for (i in 4:ncol(control_c7_freqs)){
  control_c7_freqs <- subset(control_c7_freqs, control_c7_freqs[i] >= 0 & control_c7_freqs[i] <= 1 )
}
control_c15 <-  read.table("Control//Control_wk15_haps.txt", header = TRUE)
founders <- as.vector(str_split_fixed(control_c15[1,5], ";",12))
control_c15_freqs <- cbind(control_c15[,1:3],colsplit(control_c15$adjfounderfreqs, ";", founders))
for (i in 4:ncol(control_c15_freqs)){
  control_c15_freqs <- subset(control_c15_freqs, control_c15_freqs[i] >= 0 & control_c15_freqs[i] <= 1 )
}

control_all <- rbind(control_c1_freqs, control_c7_freqs, control_c15_freqs)
control_focal <- control_all[control_all$chr == "C02" & abs(control_all$pos - general_max_pos) <= 30000, ]

moderate_c1 <-  read.table("Moderate/Mod_wk1_haps.txt", header = TRUE)
founders <- as.vector(str_split_fixed(moderate_c1[1,5], ";",12))
moderate_c1_freqs <- cbind(moderate_c1[,1:3],colsplit(moderate_c1$adjfounderfreqs, ";", founders))
for (i in 4:ncol(moderate_c1_freqs)){
  moderate_c1_freqs <- subset(moderate_c1_freqs, moderate_c1_freqs[i] >= 0 & moderate_c1_freqs[i] <= 1 )
}
moderate_c7 <-  read.table("Moderate/Mod_wk7_haps.txt", header = TRUE)
founders <- as.vector(str_split_fixed(moderate_c7[1,5], ";",12))
moderate_c7_freqs <- cbind(moderate_c7[,1:3],colsplit(moderate_c7$adjfounderfreqs, ";", founders))
for (i in 4:ncol(moderate_c7_freqs)){
  moderate_c7_freqs <- subset(moderate_c7_freqs, moderate_c7_freqs[i] >= 0 & moderate_c7_freqs[i] <= 1 )
}
moderate_c15 <-  read.table("Moderate/Mod_wk15_haps.txt", header = TRUE)
founders <- as.vector(str_split_fixed(moderate_c15[1,5], ";",12))
moderate_c15_freqs <- cbind(moderate_c15[,1:3],colsplit(moderate_c15$adjfounderfreqs, ";", founders))
for (i in 4:ncol(moderate_c15_freqs)){
  moderate_c15_freqs <- subset(moderate_c15_freqs, moderate_c15_freqs[i] >= 0 & moderate_c15_freqs[i] <= 1 )
}

moderate_all <- rbind(moderate_c1_freqs, moderate_c7_freqs, moderate_c15_freqs)
moderate_focal <- moderate_all[moderate_all$chr == "C02" & abs(moderate_all$pos - general_max_pos) <= 30000, ]

high_c1 <-  read.table("High/High_wk1_haps.txt", header = TRUE)
founders <- as.vector(str_split_fixed(high_c1[1,5], ";",12))
high_c1_freqs <- cbind(high_c1[,1:3],colsplit(high_c1$adjfounderfreqs, ";", founders))
for (i in 4:ncol(high_c1_freqs)){
  high_c1_freqs <- subset(high_c1_freqs, high_c1_freqs[i] >= 0 & high_c1_freqs[i] <= 1 )
}

high_c7 <-  read.table("High/High_wk7_haps.txt", header = TRUE)
founders <- as.vector(str_split_fixed(high_c7[1,5], ";",12))
high_c7_freqs <- cbind(high_c7[,1:3],colsplit(high_c7$adjfounderfreqs, ";", founders))
for (i in 4:ncol(high_c7_freqs)){
  high_c7_freqs <- subset(high_c7_freqs, high_c7_freqs[i] >= 0 & high_c7_freqs[i] <= 1 )
}

high_c15 <-  read.table("High/High_wk15_haps.txt", header = TRUE)
founders <- as.vector(str_split_fixed(high_c15[1,5], ";",12))
high_c15_freqs <- cbind(high_c15[,1:3],colsplit(high_c15$adjfounderfreqs, ";", founders))
for (i in 4:ncol(high_c15_freqs)){
  high_c15_freqs <- subset(high_c15_freqs, high_c15_freqs[i] >= 0 & high_c15_freqs[i] <= 1 )
}

high_all <- rbind(high_c1_freqs, high_c7_freqs, high_c15_freqs)
high_focal <- high_all[high_all$chr == "C02" & abs(high_all$pos - general_max_pos) <= 30000, ]

#colors for haps
colors <- c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231", "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE", "#008080", "#E6BEFF")




#control plot
grouped_control <- control_focal %>%
  group_by(chr, pos) %>%
  filter(n() > 1) %>%
  ungroup()

grouped_control <- grouped_control %>%
  mutate(sample_ending = sub(".*_(wk\\d+)$", "\\1", poolroot))


control_long <- grouped_control %>%
  pivot_longer(cols = starts_with("A") | starts_with("B"), names_to = "Haplotypes", values_to = "value")

control_long <- control_long %>%
  mutate(sample_ending = recode(sample_ending,
                                "wk01" = "Cycle 1",
                                "wk07" = "Cycle 7",
                                "wk15" = "Cycle 15"),sample_ending = factor(sample_ending, levels = c("Cycle 1", "Cycle 7", "Cycle 15")))

a <- ggplot(control_long, aes(x = pos, y = value, color = Haplotypes, group = interaction(poolroot,Haplotypes))) +
  geom_line() +
  facet_grid(rows = vars(sample_ending)) +
  labs(title = "A. Control Treatment",
       x = "",
       y = "Haplotype Frequency",
       color = "Haplotypes") +
  ylim(0, 1)+
  theme(legend.position = "none") +
  scale_color_manual(values = colors)

#moderate plot
grouped_moderate <- moderate_focal %>%
  group_by(chr, pos) %>%
  filter(n() > 1) %>%
  ungroup()

grouped_moderate <- grouped_moderate %>%
  mutate(sample_ending = sub(".*_(wk\\d+)$", "\\1", poolroot))


moderate_long <- grouped_moderate %>%
  pivot_longer(cols = starts_with("A") | starts_with("B"), names_to = "Haplotypes", values_to = "value")

moderate_long <- moderate_long %>%
  mutate(sample_ending = recode(sample_ending,
                                "wk01" = "Cycle 1",
                                "wk07" = "Cycle 7",
                                "wk15" = "Cycle 15"),sample_ending = factor(sample_ending, levels = c("Cycle 1", "Cycle 7", "Cycle 15")))


b <- ggplot(moderate_long, aes(x = pos, y = value, color = Haplotypes, group = interaction(poolroot,Haplotypes))) +
  geom_line() +
  facet_grid(rows = vars(sample_ending)) +
  labs(title = "B. Moderate Ethanol Treatment",
       x = "",
       y = "Haplotype Frequency",
       color = "Haplotypes") +
  ylim(0, 1)+
  theme(legend.position = "none")+
  scale_color_manual(values = colors)

#high plot
grouped_high <- high_focal %>%
  group_by(chr, pos) %>%
  filter(n() > 1) %>%
  ungroup()

grouped_high <- grouped_high %>%
  mutate(sample_ending = sub(".*_(wk\\d+)$", "\\1", poolroot))

high_long <- grouped_high %>%
  pivot_longer(cols = starts_with("A") | starts_with("B"), names_to = "Haplotypes", values_to = "value")

high_long <- high_long %>%
  mutate(sample_ending = recode(sample_ending,
                                "wk01" = "Cycle 1",
                                "wk07" = "Cycle 7",
                                "wk15" = "Cycle 15"),sample_ending = factor(sample_ending, levels = c("Cycle 1", "Cycle 7", "Cycle 15")))


c <- ggplot(high_long, aes(x = pos, y = value, color = Haplotypes, group = interaction(poolroot,Haplotypes))) +
  geom_line() +
  facet_grid(rows = vars(sample_ending)) +
  labs(title = "C. High Ethanol Treatment",
       x = "",
       y = "Haplotype Frequency",
       color = "Haplotypes") +
  ylim(0, 1) +
  theme(legend.position = "none")+
  scale_color_manual(values = colors)

#make legend

legend <- get_legend(ggplot(high_long, aes(x = pos, y = value, color = Haplotypes, group = interaction(poolroot,Haplotypes))) +
                       geom_line() +
                       facet_grid(rows = vars(sample_ending)) +
                       labs(title = "C. High Ethanol Treatment",
                            x = "Position on C11",
                            y = "Haplotype Frequency",
                            color = "Haplotypes") +
                       scale_color_manual(values = colors)+
                       ylim(0, 1)+theme(legend.position = "bottom")+
                       guides(color = guide_legend(nrow = 1)))



x_axis_label <- textGrob("Highest General Ethanol Selection Peak - chr2:108600-168600", gp = gpar(fontsize = 15))



combined_plot <- grid.arrange(
  arrangeGrob(a, b, c, ncol = 3),
  x_axis_label,legend,
  ncol = 1,
  heights = c(10, 0.1,1)  # Adjust the heights to allocate space for the legend
)

# change in B10
#control 
c_wk1 <- subset(grouped_control, sample_ending == "wk01")[,(1:15)]
c_wk15 <- subset(grouped_control, sample_ending == "wk15")[,(1:15)]
mean(c_wk15$B10) - mean(c_wk1$B10)

#moderate 
m_wk1 <- subset(grouped_moderate, sample_ending == "wk01")[,(1:15)]
m_wk15 <- subset(grouped_control, sample_ending == "wk15")[,(1:15)]
mean(m_wk15$B10) - mean(m_wk1$B10)

#high 
h_wk1 <- subset(grouped_high, sample_ending == "wk01")[,(1:15)]
h_wk15 <- subset(grouped_control, sample_ending == "wk15")[,(1:15)]
mean(h_wk15$B10) - mean(h_wk1$B10)

