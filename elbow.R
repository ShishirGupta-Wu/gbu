library(factoextra)
library(NbClust)

setwd ("C:/Users/User/Desktop/Ongoing_PROJECTS/Thapa_temp/SCR")

data <- read.table('Gene_counts_all.txt', header=T)

df <- data
df <- as.data.frame(data)

# Remove column Orthogroup
row.names(df) <- paste(df$Orthogroup, row.names(df)) 
df$Orthogroup <- NULL
head(df)



pdf(file="elbow.pdf")
fviz_nbclust(df, kmeans, method = "wss") + geom_vline(xintercept = 3, linetype = 2)+   labs(subtitle = "Elbow method")
dev.off()
