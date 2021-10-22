#Assignment3
library(DESeq2)
library(ggplot2)
library(magrittr)
library(devtools)
library(ComplexHeatmap)
library(circlize)
library(ggfortify)

#1 - R expression data and matching metadata
countdata <- read.table("./data/countdata.txt")
countdata <- countdata[,-1]

groups <- read.table("./data/groups.txt")
groups <- groups[1:90,]

#ayo what is the metadata is that the groups
#Yeah, though each application tends to need more metadata fields so make sure groups actually contains what you need

#1a
top5k <- apply(countdata, FUN=mad, MARGIN=1)
top5k <- sort(top5k, decreasing=TRUE)
top5k <- top5k[1:5000]
top5k

#1 Brian - Heirarchical
install.packages("class")
library(class)
install.packages("hclust")
install.packages("rafalib")
library(rafalib)

top5kvar <- apply(countdata, FUN=var, MARGIN=1)
top5kvar <- sort(top5kvar, decreasing=TRUE)
top5kvar <- top5kvar[1:5000]
top5kvarframe <- as.data.frame(top5kvar)

countframe = as.data.frame(countdata)
top5kcountdata <- countframe[rownames(countframe) %in% rownames(top5kvarframe),]

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x))) }

top5knorm <- as.data.frame(lapply(top5kcountdata, normalize))
colnames(top5knorm) <- groups$Combined
rownames(top5knorm) <- rownames(top5kcountdata)

t_top5knorm <- as.data.frame(t(top5knorm))

clusters <- hclust(dist(t_top5knorm))
clusterCut <- cutree(clusters, 7)

myplclust(clusters, lab.col=clusterCut, cex=0.5)

pca <- prcomp(t_top5knorm, scale=TRUE)
autoplot(pca, col=clusterCut)


hc_dend <- myplclust(clusters, lab.col=clusterCut, cex=0.5)
col_labels <- get_leaves_branches_col(hc_dend)
library(RColorBrewer)
h <- heatmap(as.matrix(clusters), Rowv=hc_dend, Colv = hc_dend, col = colorRampPalette(brewer.pal(8, "Oranges"))(12))
legend(x="topleft", legend = c("Min", "Ave", "Max"), fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))



install.packages("ggalluvial")
library(ggalluvial)

ggplot(gene_number_clusters,
       aes(y=1, axis1=Tool, axis2=Gene_Count)) +
  geom_alluvium(aes(fill=Clusters)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Survey", "Response"),
                   expand = c(0.15, 0.05))

CCPlus <- c(1, 2, 3, 2, 2, 4, 4, 3, 3, 3, 5, 3, 4, 3, 4, 3, 3, 2, 2, 4, 3, 4, 3, 3, 4, 1, 2, 5, 3, 3, 5, 1, 3, 4, 4, 5, 1, 4, 4, 4, 3, 3, 1, 6, 5, 1, 4, 6, 6, 6, 1, 4, 3, 6, 4, 3, 3, 4, 4, 6, 4, 4, 3, 4, 1, 3, 4, 2, 2, 2, 2, 4, 7, 1, 3, 3, 4, 4, 4, 2, 2, 2, 2, 3, 2, 4, 2, 3, 7, 5)
PAM <- c(1, 2, 2, 1, 1, 1, 2, 3, 2, 1, 4, 4, 2, 1, 3, 3, 2, 1, 1, 3, 1, 2, 3, 1, 1, 2, 1, 5, 3, 1, 1, 2, 1, 1, 1, 3, 4, 1, 4, 3, 3, 4, 4, 1, 1, 1, 3, 3, 3, 4, 4, 3, 4, 4, 3, 3, 1, 1, 4, 4, 4, 3, 4, 4, 3, 4, 2, 1, 2, 1, 1, 2, 4, 2, 2, 1, 1, 1, 2, 1, 2, 2, 2, 2, 1, 2, 3, 4, 1, 1)

chisq.test(table(groups$Group, groups$PAM))
chisq.test(table(groups$Experimental, groups$PAM))
chisq.test(table(groups$Month, groups$PAM))

p_values <- c(.1379, .3336, .2656, .41, .04043, .1824, .003883, .00007761, .00006631, .1806, 0, .06633)
p_adj <- p.adjust(p_values, "bonferroni")

enrichment[1,] <- c("Hclust vs treatment", 1)
enrichment[2,] <- c("CCPlus vs treatment", 1)
enrichment[3,] <- c("PAM vs treatment", 1)
enrichment[4,] <- c("Hclust vs binary treatment", 1)
enrichment[5,] <- c("CCPlus vs binary treatment", 0.48516)
enrichment[6,] <- c("PAM vs binary treatment", 1)
enrichment[7,] <- c("Hclust vs sample age", .046596)
enrichment[8,] <- c("CCPlus vs sample age", .00093132)
enrichment[9,] <- c("PAM vs sample age", .00079572)
enrichment[10,] <- c("Hclust vs CCPlus", 1)
enrichment[11,] <- c("Hclust vs PAM", 0)
enrichment[12,] <- c("CCPlus vs PAM", .79596)