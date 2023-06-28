#read in results and make S table
control <- read.table("Control_BaitER_results.txt", header = TRUE)
mod <- read.table("Moderate_BaitER_results.txt", header = TRUE)
high <- read.table("High_BaitER_results.txt", header = TRUE)
data <- data.frame(control[,1:4], mod$sigma, high$sigma)
names(data) <- c("CHR", "BP", "REF", "Control","Moderate", "High")

#read in category lists
high_only <- read.table("High_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
mod_only <- read.table("Mod_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))
con_only <- read.table("Control_only_sig_snps.txt", header =TRUE,col.names = c("CHR","BP"))

#All_Sites 
All_Sites <- data[,4:6]
#high corr
High_Specific <- na.omit(merge(data, high_only))[,4:6]
#names(high) <- c("Control_H","Moderate_H", "High_H")
#mod corr
Moderate_Specific <- na.omit(merge(data, mod_only))[,4:6]
#names(moderate) <- c("Control_M","Moderate_M", "High_M")
#control corr
Control_Specifc <- na.omit(merge(data, con_only))[,4:6]
#names(control) <- c("Control_C","Moderate_C", "High_C")

#make a table that has  corr results for each of the above + confidence intervals
tables <- list(All_Sites = All_Sites, High_Specific = High_Specific, Moderate_Specific = Moderate_Specific, Control_Specifc = Control_Specifc)

result_table <- data.frame(
  Category = character(),
  Group1 = character(),
  Group2 = character(),
  Correlation = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
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
        CI_Upper = ci[2]
      ))
    }
  }
}

result_table$Category <- gsub("_", " ", result_table$Categor)

write.table (result_table,file="Pearson_corr_results.txt",quote=FALSE,row.names = FALSE, col.names = TRUE, sep = "\t")
