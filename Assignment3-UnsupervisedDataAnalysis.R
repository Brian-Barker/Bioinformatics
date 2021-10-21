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