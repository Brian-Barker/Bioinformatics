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

#1 Brian - Gaussian Mixture Models
install.packages("mclust")
library(mclust)
top5kvar <- apply(countdata, FUN=var, MARGIN=1)
top5kvar <- sort(top5kvar, decreasing=TRUE)
top5kvar <- top5kvar[1:5000]
top5kvar
top5kvarframe = as.data.frame(top5kvar)
mod <- mclust::Mclust(top5kvar)
dr <- MclustDR(mod)

clPairs(top5kvar, groups$Experimental)

#Paige - PAM Clustering
countframe = as.data.frame(countdata)
top5kcountdata <- countframe[rownames(countframe) %in% rownames(top5kvarframe),]

install.packages(c("cluster", "factoextra"))
library(cluster)
library(factoextra)

pamx <- pam(top5kcountdata, 5, metric = "euclidean")
plot(pamx)

top10 <- apply(countdata, FUN=mad, MARGIN=1)
top10 <- sort(top1k, decreasing=TRUE)
top10 <- top10[1:10]
top10varframe = as.data.frame((top10var))
top10countdata <- countframe[rownames(countframe) %in% rownames(top10varframe),]
top10norm <- as.data.frame(lapply(top10countdata, normalize))
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x))) }
colnames(top10norm) <- groups$Combined
rownames(top10norm) <- rownames(top10countdata)

t_top10norm <- as.data.frame(t(top10norm))



top100 <- apply(countdata, FUN=mad, MARGIN=1)
top100 <- sort(top1k, decreasing=TRUE)
top100 <- top100[1:100]
top100varframe = as.data.frame((top100var))
top100countdata <- countframe[rownames(countframe) %in% rownames(top100varframe),]
top100norm <- as.data.frame(lapply(top100countdata, normalize))
colnames(top100norm) <- groups$Combined
rownames(top100norm) <- rownames(top100countdata)

t_top100norm <- as.data.frame(t(top100norm))

top5k <- apply(countdata, FUN=mad, MARGIN=1)
top5k <- sort(top5k, decreasing=TRUE)
top5k <- top5k[1:5000]
top5kvarframe = as.data.frame((top5kvar))
top5kcountdata <- countframe[rownames(countframe) %in% rownames(top5kvarframe),]
top5knorm <- as.data.frame(lapply(top5kcountdata, normalize))
colnames(top5knorm) <- groups$Combined
rownames(top5knorm) <- rownames(top5kcountdata)

t_top5knorm <- as.data.frame(t(top5knorm))


top10k <- apply(countdata, FUN=mad, MARGIN=1)
top10k <- sort(top10k, decreasing=TRUE)
top10k <- top10k[1:10000]
top10kvarframe = as.data.frame((top10kvar))
top10kcountdata <- countframe[rownames(countframe) %in% rownames(top10kvarframe),]
top10knorm <- as.data.frame(lapply(top10kcountdata, normalize))
colnames(top10knorm) <- groups$Combined
rownames(top10knorm) <- rownames(top10kcountdata)

t_top10knorm <- as.data.frame(t(top10knorm))


top1k <- apply(countdata, FUN=mad, MARGIN=1)
top1k <- sort(top1k, decreasing=TRUE)
top1k <- top1k[1:1000]
top1kvarframe = as.data.frame((top1kvar))
top1kcountdata <- countframe[rownames(countframe) %in% rownames(top1kvarframe),]
top1knorm <- as.data.frame(lapply(top1kcountdata, normalize))
colnames(top1knorm) <- groups$Combined
rownames(top1knorm) <- rownames(top1kcountdata)

t_top1knorm <- as.data.frame(t(top1knorm))

plotting <- t_top5knorm
#5k   2
#100  4
#1k   2
#10k  2
#10   4

pam.res <- pam(plotting, k=4)
fviz_cluster(pam.res, geom = "point", ellipse.type = "norm")

fviz_nbclust(plotting, pam, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

fviz_nbclust(plotting, pam, method = "silhouette")+
  labs(subtitle = "Silhouette method")


library(RColorBrewer)
h <- heatmap(as.matrix(pam.res), col = colorRampPalette(brewer.pal(8, "Oranges"))(12))
legend(x="topleft", legend = c("Min", "Ave", "Max"), fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))


heatmap(as.matrix(top5knorm, name = "mat", row_split = paste0("pam", pam.res$clustering))

mat <- as.matrix(top5knorm), nrow = 20)
PAM.hm(as.matrix(top5knorm), cluster.number = 3)

       