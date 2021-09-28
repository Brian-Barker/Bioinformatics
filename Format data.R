countdata <- subset(data, select=-c(Chr, Start, End, strand, Length))
countdata <- countdata[,-1]
rownames(countdata) <- countdata[,1]

librarySizes <- colSums(countdata)
barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, 
        main="Barplot of library sizes")