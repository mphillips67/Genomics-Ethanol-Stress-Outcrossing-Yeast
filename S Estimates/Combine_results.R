#combine S estimates into one file

#read in results
control <- read.table("Control_BaitER_results.txt", header = TRUE)
mod <- read.table("Moderate_BaitER_results.txt", header = TRUE)
high <- read.table("High_BairER_results.txt", header = TRUE)

#combine
control_s <- control$sigma
mod_s <- mod$sigma
high_s <- high$sigma

data <- cbind(control[,1:3],control_s,mod_s,high_s)

write.table (data,file="BaitER_results_all.txt",quote=FALSE,row.names = FALSE, col.names = TRUE, sep = "\t")
