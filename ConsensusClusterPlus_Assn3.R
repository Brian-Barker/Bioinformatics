##########################
#
# ConsensusClusterPlus
#
##########################

#install packages
library(DESeq2)
library(ggplot2)
library(magrittr)
library(devtools)
library(ComplexHeatmap)
library(circlize)
library(ggfortify)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ConsensusClusterPlus")

#load data
countdata <- read.table("./data/countdata.txt")
countdata <- countdata[,1:90]
cols_ <- colnames(countdata)
rows_ <- rownames(countdata)
countdata <- matrix(unlist(countdata), ncol=90, nrow=lengths(countdata))

groups <- read.table("./data/groups_nospace.txt")
groups <- groups[1:90,]

cols_actual <- groups[1:90, 3]
rownames(countdata) = rows_
colnames(countdata) = cols_

#format specifically for 5k
top_genes <- apply(countdata, FUN=mad, MARGIN=1)
countdata_5k = countdata[rev(order(top_genes))[1:5000],]

#normalize the data
countdata_5k = sweep(countdata_5k, 1, apply(countdata_5k, 1,median, na.rm=T))

library(ConsensusClusterPlus)
title=tempdir()
results_5k = ConsensusClusterPlus(countdata_5k,maxK=10,reps=1000,pItem=0.8,pFeature=1,title="10_5k",clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png", writeTable=TRUE)
tmp=read.csv("./10_5k/10_5k.k=6.consensusClass.csv", header=FALSE)
tmp=cbind(cols_actual, tmp)
write.csv(tmp,"./10_5k/10_5k_actual.k=6.consensusClass.csv")

#format specifically for 10 genes
countdata_10 = countdata[rev(order(top_genes))[1:10],]
#normalize the data
countdata_10 = sweep(countdata_10, 1, apply(countdata_10, 1,median, na.rm=T))
#run consensusclusterplus
title=tempdir()
results_10 = ConsensusClusterPlus(countdata_10,maxK=10,reps=1000,pItem=0.8,pFeature=1,title="10_10",clusterAlg="hc",distance="pearson",seed=1262118388.71278,plot="png", writeTable=TRUE)

#format specifically for 100 genes
countdata_100 = countdata[rev(order(top_genes))[1:100],]
#normalize the data
countdata_100 = sweep(countdata_100, 1, apply(countdata_100, 1,median, na.rm=T))
#run consensusclusterplus
title=tempdir()
results_100 = ConsensusClusterPlus(countdata_100,maxK=10,reps=1000,pItem=0.8,pFeature=1,title="10_100",clusterAlg="hc",distance="pearson",seed=1262118388.71277,plot="png", writeTable=TRUE)

#format specifically for 1000 genes
countdata_1000 = countdata[rev(order(top_genes))[1:1000],]
#normalize the data
countdata_1000 = sweep(countdata_1000, 1, apply(countdata_1000, 1,median, na.rm=T))
#run consensusclusterplus
title=tempdir()
results_1000 = ConsensusClusterPlus(countdata_1000,maxK=10,reps=1000,pItem=0.8,pFeature=1,title="10_1000",clusterAlg="hc",distance="pearson",seed=1262118388.71276,plot="png", writeTable=TRUE)

#format specifically for 10000 genes
countdata_10000 = countdata[rev(order(top_genes))[1:10000],]
#normalize the data
countdata_10000 = sweep(countdata_10000, 1, apply(countdata_10000, 1,median, na.rm=T))
#run consensusclusterplus
title=tempdir()
results_10000 = ConsensusClusterPlus(countdata_10000,maxK=10,reps=1000,pItem=0.8,pFeature=1,title="10_10000",clusterAlg="hc",distance="pearson",seed=1262118388.71275,plot="png", writeTable=TRUE)

#Heatmap
library(RColorBrewer)
consensus_ <- results_5k[[7]][["consensusMatrix"]]
colnames(consensus_) = cols_actual
rownames(consensus_) = cols_actual

heatmap(consensus_, col = colorRampPalette(brewer.pal(8, "Blues"))(12))
legend(x="topleft", legend = c("Min", "Ave", "Max"), fill=colorRampPalette(brewer.pal(8, "Blues"))(3))

#STATS BAYBAY
tmp=read.csv("./10_5k/10_5k.k=7.consensusClass.csv", header=FALSE)
df <- tmp[,2]
