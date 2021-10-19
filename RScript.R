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

library(M3C)
tsne(logcounts, labels=groups[,4], text=groups[,2])


require(DOSE)
require(clusterProfiler)
orange <- search_kegg_organism("cic", by = "kegg_code")
dim(orange)
head(orange)
data(geneList, package="DOSE")
k <- enrichMKEGG(gene = "CICLE_v10004626mg",
               organism = 'cic',
               pvalueCutoff = 0.05)
head(k)

## this didn't upload, but 1 and 5 i think
data_file <- "./data/GSE167171_AllCounts.txt"
d <- read.table(data_file, TRUE, sep="")

#####################
# 1
#####################
#1c
columns <- ncol(d) - 7
rows <-nrow(d)
cat("There are a total of ", toString(columns), " samples in the data.")
cat("Additionally, there are a total of ", toString(rows), "genes that were differentially expressed in the sample")
#1c - what is the variation in the data
v <- d[, 7:ncol(d)]
vary <- var(v)
range(vary)
#1c - per gene expression ranges
rowranges <- numeric()
for (i in 1:nrow(d)) {
        rowranges <- c(rowranges, diff(range(d[i, 7:ncol(d)])))
}
ct <- fivenum(rowranges)
upper <- ct[4] +(ct[4]-ct[2])*1.5
lower <- ct[2] -(ct[4]-ct[2])*1.5
rowranges <- rowranges[rowranges<=upper]
dens <- density(rowranges)
autoplot(dens, xlab="Range Values of Gene Expression", ylab="Density", main="Per-Gene Expression Ranges")

#################
# 4
#################

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = groups, design = ~Group)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = groups, design = ~Group + Month)
dds <- DESeq(dds)
results_dds <- results(dds)
results_dds <- lfcShrink( dds, coef=2, res = results_dds)

dds.df <- results_dds %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Gene") %>%
        dplyr::mutate(threshold = padj < 0.05) %>%
        dplyr::arrange(dplyr::desc(log2FoldChange))

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)

#building heatmap
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "blue"))
col_fun(seq(-3, 3))

normz = results_dds$padj <= 0.01 & results_dds$baseMean >= 0.05 & 
        abs(results_dds$log2FoldChange) >= 0; normz[is.na(normz)] = FALSE
dds_m = counts(dds, normalized = TRUE)
dds_m = dds_m[normz, ]

heat = Heatmap(t(scale(t(dds_m))), name = "z-score",
               top_annotation = HeatmapAnnotation(
                       TreatmentGroups = colData(dds)$Group,
                       sizeFactor = anno_points(colData(dds)$sizeFactor)
               ),
               show_row_names = FALSE, show_column_names = FALSE,
               column_title = paste0(sum(normz), " Significant Genes Assuming FDR = 0.05"),
               show_row_dend = FALSE)+
        Heatmap(results_dds$log2FoldChange[normz], show_row_names = FALSE, width = unit(7, "mm"),
                name = "log2FoldChange", show_column_names = FALSE, col=col_fun)
heat = draw(heat, merge_legend = TRUE)
heat