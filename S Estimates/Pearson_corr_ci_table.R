#make table with 95 CI and SNP counts for correlations between top candidates for each category. requires cmh files, s estimate files, and category files. 
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

#All_Sites 
All_Sites <- data[,4:6]

#general Lab
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
General_Lab_Selection <- na.omit(merge(data, general_lab_msm_range))[,4:6]

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
General_Ethanol_Selection <- na.omit(merge(data, general_eth_msm_range))[,4:6]

#high corr
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
High_Specific <- na.omit(merge(data, high_eth_msm_range))[,4:6]

#mod corr
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
Moderate_Specific <- na.omit(merge(data, mod_eth_msm_range))[,4:6]

#control corr
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
Control_Specific <- na.omit(merge(data, control_eth_msm_range))[,4:6]


#make a table that has  corr results for each of the above + confidence intervals
tables <- list(All_Sites = All_Sites,General_Lab_Selection= General_Lab_Selection, General_Ethanol_Selection =General_Ethanol_Selection , High_Specific = High_Specific, Moderate_Specific = Moderate_Specific, Control_Specific = Control_Specific)

result_table <- data.frame(
  Category = character(),
  Group1 = character(),
  Group2 = character(),
  Correlation = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  SNP_count = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each table
for (table_name in names(tables)) {
  # Get the current table
  current_table <- tables[[table_name]]
  
  # Calculate pairwise Pearson correlations
  cor_matrix <- cor(current_table, use = "complete.obs", method = "pearson")
  
  # Loop through All_Sites pairwise combinations of columns
  for (i in 1:(ncol(current_table) - 1)) {
    for (j in (i + 1):ncol(current_table)) {
      # Get column names
      var1 <- colnames(current_table)[i]
      var2 <- colnames(current_table)[j]
      
      # Calculate correlation and confidence interval
      correlation <- cor_matrix[i, j]
      ci <- cor.test(current_table[[var1]], current_table[[var2]])$conf.int
      
      # Add results to the table
      result_table <- rbind(result_table, data.frame(
        Category = table_name,
        Group1 = var1,
        Group2 = var2,
        Correlation = correlation,
        CI_Lower = ci[1],
        CI_Upper = ci[2],
        SNP_count = nrow(current_table)
      ))
    }
  }
}

result_table$Category <- gsub("_", " ", result_table$Categor)

write.table (result_table,file="Pearson_corr_results.txt",quote=FALSE,row.names = FALSE, col.names = TRUE, sep = "\t")
