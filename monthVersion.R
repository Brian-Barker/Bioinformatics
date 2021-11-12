month2groups <- matrix(nrow = length(month2counts), ncol = 3)
month4groups <- matrix(nrow = length(month4counts), ncol = 3)
month8groups <- matrix(nrow = length(month8counts), ncol = 3)
month2groups[,1] <- substr(colnames(month2counts),4,4)
month4groups[,1] <- substr(colnames(month4counts),4,4)
month8groups[,1] <- substr(colnames(month8counts),4,4)
month2groups[,2] <- substr(colnames(month2counts),5,5)
month4groups[,2] <- substr(colnames(month4counts),5,5)
month8groups[,2] <- substr(colnames(month8counts),5,5)
month2groups[,3] <- ifelse (substr(colnames(month2counts),4,4) == "B", "Control", "Experimental")
month4groups[,3] <- ifelse (substr(colnames(month4counts),4,4) == "B", "Control", "Experimental")
month8groups[,3] <- ifelse (substr(colnames(month8counts),4,4) == "B", "Control", "Experimental")

colnames(month2groups) <- c("Group", "ID", "Experimental")
colnames(month4groups) <- c("Group", "ID", "Experimental")
colnames(month8groups) <- c("Group", "ID", "Experimental")

month2counts <- countdata[substr(colnames(countdata),3,3) == 2]
month4counts <- countdata[substr(colnames(countdata),3,3) == 4]
month8counts <- countdata[substr(colnames(countdata),3,3) == 8]

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x))) }

# --- PCA ---
library(ggplot2)
library(ggfortify)
create_pca <- function(counts_, groups_, title_) {
  counts_[counts_ == 0] <- 1
  logcounts <- log2(counts_)
  pcDat <- prcomp(t(logcounts))
  autoplot(pcDat, data = groups_, colour="Group", shape="Experimental") + ggtitle(title_)
}

create_pca(month2counts, month2groups, "Month 2")
create_pca(month4counts, month4groups, "Month 4")
create_pca(month8counts, month8groups, "Month 8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
library(DESeq2)
library(EnhancedVolcano)

create_volcano <- function(counts_, groups_, title_) {
  dds <- DESeqDataSetFromMatrix(countData = counts_, colData = groups_, design = ~Experimental)
  dds <- DESeq(dds, betaPrior = FALSE)
  
  
  res <- results(dds,
                 contrast = c('Experimental','Experimental','Control'))
  res <- lfcShrink(dds,
                   contrast = c('Experimental','Experimental','Control'), res=res, type = 'normal')
  EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  title = title_,
                  pCutoff = 10e-4,
                  FCcutoff = 0.3,
                  xlim = c(-1, 1),
                  ylim = c(0, 8))
}

create_volcano(month2counts, month2groups, "Month 2")
create_volcano(month4counts, month4groups, "Month 4")
create_volcano(month8counts, month8groups, "Month 8")