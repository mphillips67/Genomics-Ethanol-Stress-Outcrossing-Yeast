library(gridExtra)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(grid)
library(ggcorrplot)

#read in results and make S table
control <- read.table("Control_BaitER_results.txt", header = TRUE)
mod <- read.table("Moderate_BaitER_results.txt", header = TRUE)
high <- read.table("High_BaitER_results.txt", header = TRUE)
data <- data.frame(control[,1:4], mod$sigma, high$sigma)
names(data) <- c("CHR", "BP", "REF", "C-populations","M-populations", "H-populations")

#read in category lists
high_only <- read.table("High_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
mod_only <- read.table("Mod_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
control_only <- read.table("Control_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
general_lab  <- read.table("Sig_SNPs_Outcrossing.txt", header =TRUE,col.names = c("CHR","BP"))
general_eth <- read.table("Sig_SNPs_Shared_High_Mod_Overlap.txt", header =TRUE,col.names = c("CHR","BP"))

#genome wide s r2 plot
corr_all <- cor(data[,4:6], use = "complete.obs", method = "pearson") #cor matrix
all <- data[,4:6]
result_table_all <- data.frame(
  Variable1 = character(),
  Variable2 = character(),
  Correlation = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  stringsAsFactors = FALSE
) #I want a table with the confidence intervals associated with cor matrix

for (i in 1:(ncol(all) - 1)) {
  for (j in (i + 1):ncol(all)) {
    # Get column names
    var1 <- colnames(all)[i]
    var2 <- colnames(all)[j]
    
    # Calculate correlation and confidence interval
    correlation <- corr_all[i, j]
    ci <- cor.test(all[[var1]], all[[var2]])$conf.int
    
    # Add results to the table
    result_table_all <- rbind(result_table_all, data.frame(
      Variable1 = var1,
      Variable2 = var2,
      Correlation = correlation,
      CI_Lower = ci[1],
      CI_Upper = ci[2]
    ))
  }
}

CI_all <- paste(round(result_table_all$CI_Lower,2), round(result_table_all$CI_Upper,2), sep = "-")

a <- ggcorrplot(corr_all, digits= 2, method = "square", colors = c("#E46726", "white", "#6D9EC1"), outline.col = "white", lab = TRUE, type = "lower", title = "A. All Polymorphic Sites") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_text(aes(label=CI_all), vjust = 3,size=3)

#general lab s r2 plot
lab <- na.omit(merge(data, general_lab))[,4:6]
corr_lab <- cor(lab, use = "complete.obs", method = "pearson")
result_table_lab <- data.frame(
  Variable1 = character(),
  Variable2 = character(),
  Correlation = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  stringsAsFactors = FALSE
) #I want a table with the confidence intervals associated with cor matrix

for (i in 1:(ncol(lab) - 1)) {
  for (j in (i + 1):ncol(lab)) {
    # Get column names
    var1 <- colnames(lab)[i]
    var2 <- colnames(lab)[j]
    
    # Calculate correlation and confidence interval
    correlation <- corr_lab[i, j]
    ci <- cor.test(lab[[var1]], lab[[var2]])$conf.int
    
    # Add results to the table
    result_table_lab <- rbind(result_table_lab, data.frame(
      Variable1 = var1,
      Variable2 = var2,
      Correlation = correlation,
      CI_Lower = ci[1],
      CI_Upper = ci[2]
    ))
  }
}

CI_lab <- paste(round(result_table_lab$CI_Lower,2), round(result_table_lab$CI_Upper,2), sep = "-")

b <- ggcorrplot(corr_lab, digits= 2,method = "square", colors = c("#E46726", "white", "#6D9EC1"), outline.col = "white", lab = TRUE, type = "lower",title = "B. General Lab Selection Candidates")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_text(aes(label=CI_lab), vjust = 3,size=3)

#general ethanol s r2 plot
ethanol <- na.omit(merge(data, general_eth))[,4:6]
corr_ethanol <- cor(ethanol, use = "complete.obs", method = "pearson")
result_table_ethanol <- data.frame(
  Variable1 = character(),
  Variable2 = character(),
  Correlation = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  stringsAsFactors = FALSE
) #I want a table with the confidence intervals associated with cor matrix

for (i in 1:(ncol(ethanol) - 1)) {
  for (j in (i + 1):ncol(ethanol)) {
    # Get column names
    var1 <- colnames(ethanol)[i]
    var2 <- colnames(ethanol)[j]
    
    # Calculate correlation and confidence interval
    correlation <- corr_ethanol[i, j]
    ci <- cor.test(ethanol[[var1]], ethanol[[var2]])$conf.int
    
    # Add results to the table
    result_table_ethanol <- rbind(result_table_ethanol, data.frame(
      Variable1 = var1,
      Variable2 = var2,
      Correlation = correlation,
      CI_Lower = ci[1],
      CI_Upper = ci[2]
    ))
  }
}

CI_ethanol <- paste(round(result_table_ethanol$CI_Lower,2), round(result_table_ethanol$CI_Upper,2), sep = "-")

c <- ggcorrplot(corr_ethanol, digits= 2,method = "square", colors = c("#E46726", "white", "#6D9EC1"), outline.col = "white", lab = TRUE, type = "lower",title = "C. General Ethanol Selection Candidates")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_text(aes(label=CI_ethanol), vjust = 3,size=3)




#high only s r2 plot 
high <- na.omit(merge(data, high_only))[,4:6]
corr_high <- cor(high, use = "complete.obs", method = "pearson")
result_table_high <- data.frame(
  Variable1 = character(),
  Variable2 = character(),
  Correlation = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  stringsAsFactors = FALSE
) #I want a table with the confidence intervals associated with cor matrix

for (i in 1:(ncol(high) - 1)) {
  for (j in (i + 1):ncol(high)) {
    # Get column names
    var1 <- colnames(high)[i]
    var2 <- colnames(high)[j]
    
    # Calculate correlation and confidence interval
    correlation <- corr_high[i, j]
    ci <- cor.test(high[[var1]], high[[var2]])$conf.int
    
    # Add results to the table
    result_table_high <- rbind(result_table_high, data.frame(
      Variable1 = var1,
      Variable2 = var2,
      Correlation = correlation,
      CI_Lower = ci[1],
      CI_Upper = ci[2]
    ))
  }
}

CI_high <- paste(round(result_table_high$CI_Lower,2), round(result_table_high$CI_Upper,2), sep = "-")

d <- ggcorrplot(corr_high, digits= 2,method = "square", colors = c("#E46726", "white", "#6D9EC1"), outline.col = "white", lab = TRUE, type = "lower",title = "D. High Ethanol Specific Candidates")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_text(aes(label=CI_high), vjust = 3,size=3)

#mod only s r2 plot
moderate <- na.omit(merge(data, mod_only))[,4:6]
corr_mod <- cor(moderate, use = "complete.obs", method = "pearson")
result_table_mod <- data.frame(
  Variable1 = character(),
  Variable2 = character(),
  Correlation = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  stringsAsFactors = FALSE
) #I want a table with the confidence intervals associated with cor matrix

for (i in 1:(ncol(moderate) - 1)) {
  for (j in (i + 1):ncol(moderate)) {
    # Get column names
    var1 <- colnames(moderate)[i]
    var2 <- colnames(moderate)[j]
    
    # Calculate correlation and confidence interval
    correlation <- corr_mod[i, j]
    ci <- cor.test(moderate[[var1]], moderate[[var2]])$conf.int
    
    # Add results to the table
    result_table_mod <- rbind(result_table_mod, data.frame(
      Variable1 = var1,
      Variable2 = var2,
      Correlation = correlation,
      CI_Lower = ci[1],
      CI_Upper = ci[2]
    ))
  }
}
CI_mod <- paste(round(result_table_mod$CI_Lower,2), round(result_table_mod$CI_Upper,2), sep = "-")

e <- ggcorrplot(corr_mod, digits= 2,method = "square",colors = c("#E46726", "white", "#6D9EC1"), outline.col = "white", lab = TRUE, type = "lower", title = "E. Moderate Ethanol Specific Candidates") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_text(aes(label=CI_mod), vjust = 3,size=3)

#control only s r2 plot
control <- na.omit(merge(data, control_only))[,4:6]
corr_con <- cor(control, use = "complete.obs", method = "pearson")
result_table_con <- data.frame(
  Variable1 = character(),
  Variable2 = character(),
  Correlation = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  stringsAsFactors = FALSE
) #I want a table with the confidence intervals associated with cor matrix

for (i in 1:(ncol(control) - 1)) {
  for (j in (i + 1):ncol(control)) {
    # Get column names
    var1 <- colnames(control)[i]
    var2 <- colnames(control)[j]
    
    # Calculate correlation and confidence interval
    correlation <- corr_con[i, j]
    ci <- cor.test(control[[var1]], control[[var2]])$conf.int
    
    # Add results to the table
    result_table_con <- rbind(result_table_con, data.frame(
      Variable1 = var1,
      Variable2 = var2,
      Correlation = correlation,
      CI_Lower = ci[1],
      CI_Upper = ci[2]
    ))
  }
}
CI_con <- paste(round(result_table_con$CI_Lower,2), round(result_table_con$CI_Upper,2), sep = "-")
f <- ggcorrplot(corr_con, digits= 2,method = "square",colors = c("#E46726", "white", "#6D9EC1"), outline.col = "white", lab = TRUE, type = "lower",
                title = "F. Control Specific Candidates") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_text(aes(label=CI_con), vjust = 3,size=3)


grid.arrange(a,b,c,d,e,f,nrow=2)


