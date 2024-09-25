install.packages("e1071")


#
#set.seed(12345)
library(e1071)

inputFile="diffGeneExp.txt"     #

source("geoFRG12.msvmRFE.R")

#
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.data.frame(t(data))
#
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=cbind(group, data)
data$group=factor(data$group, levels=c("Control","Treat"))

#RFE(data, k=10, halve.above=50)
nfold=5
geneNum=nrow(data)
folds=rep(1:nfold, len=geneNum)[sample(geneNum)]
folds=lapply(1:nfold, function(x) which(folds == x))
results=lapply(folds, svmRFE.wrap, data, k=10, halve.above=50)

#????.features=WriteFeatures(results, data, save=F)
#????te.table(top.features, file="feature_svm.txt", sep="\t", quote=F,row.names=F)

#????tsweep=lapply(1:5, FeatSweep.wrap, results, data)

#??È¡info=min(prop.table(table(data[,1])))
errors=sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
 
#???Æ(file="errors.pdf", width=5, height=5)
PlotErrors(errors, no.info=no.info)
dev.off()

#???Æ(file="accuracy.pdf", width=5, height=5)
Plotaccuracy(1-errors, no.info=no.info)
dev.off()

#????tureGenes=top.features[1:which.min(errors),1,drop=F]
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)


##