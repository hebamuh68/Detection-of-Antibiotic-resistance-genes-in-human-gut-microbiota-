# Visualization - Data plots for selected GEO samples
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE54070", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL18168", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)

# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

#====================================box-and-whisker plot============================================
png(filename="whole_data_plots/box-and-whisker.png",width=850)
par(mar=c(7,4,2,1))
title <- paste ("GSE54070", "/", annotation(gset), sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
dev.off()

#====================================expression value distribution plot============================================
png(filename="whole_data_plots/expression value.png",width=850)
par(mar=c(4,4,2,1))
title <- paste ("GSE54070", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)
dev.off()

#=================================== mean-variance trend =============================
png(filename="whole_data_plots/mean-variance.png",width=850)
ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE54070")
dev.off()