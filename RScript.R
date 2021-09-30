data_file <- "data/GSE167171_AllCounts.txt"

data <- read.table(data_file, TRUE, sep="")

countdata <- read.table("./data/countdata.txt")
countdata <- countdata[,-1]
rownames(countdata) <- countdata[,1]

groups <- read.table("./data/groups.txt")
groups <- groups[1:90,]

dim(countdata)
keep <- rowSums(countdata) > 5
countdata <- countdata[keep,]
dim(countdata)

librarySizes <- colSums(countdata)
barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, 
        main="Barplot of library sizes")

logcounts <- log2(countdata + 1)
boxplot(logcounts,
        xlab="",
        ylab="Log2(Counts)",
        las=2)
abline(h=median(as.matrix(logcounts)), col = "blue")


library(ggplot2)
library(ggfortify)
rlogcounts <- log2(counts)
pcDat <- prcomp(t(logcounts))
autoplot(pcDat, data = groups, colour="Group")
