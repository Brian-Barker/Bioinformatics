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
autoplot(pcDat, data = groups, colour="Month", alpha="Experimental", shape="Experimental") + scale_alpha_discrete(range=c(1,.5))

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = groups, design = ~Experimental + Month)
dds <- DESeq(dds, betaPrior = FALSE)


res <- results(dds,
               contrast = c('Experimental','Experimental','Control'))
res <- lfcShrink(dds,
                 contrast = c('Experimental','Experimental','Control'), res=res, type = 'normal')
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 10e-4,
                FCcutoff = 0.3,
                xlim = c(-1, 1),
                ylim = c(0, 8))

differentlyExpressedGenes <- res[which(res[["pvalue"]] < 10e-3 & res[["log2FoldChange"]] > .3),]
write.table(differentlyExpressedGenes, "./data/DifferentlyExpressedGenes.txt")