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

#1a
top5k <- apply(countdata, FUN=var, MARGIN=1)
top5k <- sort(top5k, decreasing=TRUE)
top5k <- top5k[1:5000]
