install.packages("ggplot2")

data_file <- "GSE167171_AllCounts.txt/GSE167171_AllCounts.txt"

d <- read.table(data_file, TRUE, sep="")

#control + enter to assign objects
print(d)

library(ggplot2)
ggplot(data=d, aes(x="MI2C5", y="Length"))
