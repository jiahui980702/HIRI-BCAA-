if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")


#
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

expFile="normalize_GSE23649.txt"     
gene="SLC7A5"                
gmtFile="c2.cp.kegg.symbols.gmt"  


#
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    nL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#????FC=sort(logFC, decreasing=T)
genes=names(logFC)

#??ȡ=read.gmt(gmtFile)

#GSEA????GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_SLC7A5_KEGG.txt",sep="\t",quote=F,row.names = F)
	
#????mNum=6     #չʾnrow(kkTab)>=termNum){
	showTerm=row.names(kkTab)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title=gene)
	pdf(file="GSEA_SLC7A5_KEGG.pdf", width=7.5, height=5.5)
	print(gseaplot)
	dev.off()
}









##==================================================
######VideFile="normalize_GSE23649.txt"     #??
gene="SLC7A5"               File="c5.go.symbols.gmt"     #??w?
#read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#ȥ??up=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#????ata[,data[gene,]<median(data[gene,]),drop=F]     #?ͱ???
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]  
meanL=rowMeans(dataL)owMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#????logFort(logFC, decreasing=T)
genes=names(logFC)

#??ȡ???.gmt(gmtFile)

#GSEA??
#GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_SLC7A5_GO.txt",sep="\t",quote=F,row.names = F)

#????GSEA??ͼ=6     #չʾͨ�(kkTab)>=termNum){
  showTerm=row.names(kkTab)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title=gene)
  pdf(file="GSEA_SLC7A5_GO.pdf", width=7.5, height=5.5)
  print(gseaplot)
  dev.off()
}


