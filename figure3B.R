install.packages("colorspace")
install.packages("stringi")
install.packages("ggplot2")
install.packages("circlize")
install.packages("RColorBrewer")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("ComplexHeatmap")


#
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05     #pֵ????????
qvalueFilter=1        #????????pֵ????????

#
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}


rt=read.table("diffGeneExp.txt", header=T, sep="\t", check.names=F)     #

#
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
#
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#
pdf(file="barplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
dev.off()

#
pdf(file="bubble.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
dev.off()

