

requiredPackages = c('tidyverse')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

library(tidyverse)

data <- read.table('Final_removed_cluster_counts_for_plot.txt', header=T)
df <- data
df <- as.data.frame(data)
df


pdf(file="cluster_loss_chart.pdf")
ggplot(data = df, aes(x = Threshold, y = Deleted_clusters)) +   geom_line(color = 'navy', size = 1) + scale_x_continuous(breaks=seq(0, 750, by = 50))
dev.off()


