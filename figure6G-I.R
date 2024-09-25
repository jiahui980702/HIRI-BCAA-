if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("GSEABase")
BiocManager::install("GSVA")

#install.packages("ggpubr")

rm(list = ls())   
#
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase) 
library(GSVA)

gene="SLC7A5"      
expFile="normalize_GSE23649.txt"              
gmtFile="c2.cp.kegg.symbols.gmt"     


#
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#GSVA
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
ssgseaScore=ssgseaScore[order(apply(ssgseaScore,1,sd),decreasing=T),]
ssgseaScore=ssgseaScore[1:50,]

#Name=colnames(data)[data[gene,]<median(data[gene,])]       #?ͱ?hName=colnames(data)[data[gene,]>=median(data[gene,])]     #?߱?Score=ssgseaScore[,lowName]
highScore=ssgseaScore[,highName]
data=cbind(lowScore, highScore)
conNum=ncol(lowScore)
treatNum=ncol(highScore)
Type=c(rep("Control",conNum), rep("Treat",treatNum))

#????Tab=data.frame()
for(i in row.names(data)){
	test=t.test(data[i,] ~ Type)
	pvalue=test$p.value
	t=test$statistic
	Sig=ifelse(pvalue>0.05, "Not", ifelse(t>0,"Up","Down"))
	outTab=rbind(outTab, cbind(Pathway=i, t, pvalue, Sig))
}

#????(file="barplot_SLC7A5_KEGG.pdf", width=12, height=9)
outTab$t=as.numeric(outTab$t)
outTab$Sig=factor(outTab$Sig, levels=c("Down", "Not", "Up"))
gg1=ggbarplot(outTab, x="Pathway", y="t", fill = "Sig", color = "white",
		palette=c("green","black","#E41A1C"), sort.val = "asc", sort.by.groups = T,
		rotate=TRUE, legend="right", title=gene,
		xlab="Term", ylab="t value of GSVA score",  legend.title="Group", x.text.angle=60)
print(gg1)
dev.off()




#========================================
######Vide
gene="SLC7A5"     
expFile="normalize_GSE23649.txt"            
gmtFile="c5.go.symbols.gmt"   w?r#
ead.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#ȥ
#up=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#??ȡeSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#GS
#GSVAseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#?Դ?malize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
ssgseaScore=ssgseaScore[order(apply(ssgseaScore,1,sd),decreasing=T),]
ssgseaScore=ssgseaScore[1:30,]

#????=colnames(data)[data[gene,]<median(data[gene,])]       #?ͱ?????????Ʒ
highName=colnames(data)[data[gene,]>=median(data[gene,])]     #?߱?????????Ʒ
lowScore=ssgseaScore[,lowName]
highScore=ssgseaScore[,highName]
data=cbind(lowScore, highScore)
conNum=ncol(lowScore)
treatNum=ncol(highScore)
Type=c(rep("Control",conNum), rep("Treat",treatNum))

#????????data.frame()
for(i in row.names(data)){
  test=t.test(data[i,] ~ Type)
  pvalue=test$p.value
  t=test$statistic
  Sig=ifelse(pvalue>0.05, "Not", ifelse(t>0,"Up","Down"))
  outTab=rbind(outTab, cbind(Pathway=i, t, pvalue, Sig))
}

#??????״e="barplot_GO_SLC7A5.pdf", width=12, height=9)
outTab$t=as.numeric(outTab$t)
outTab$Sig=factor(outTab$Sig, levels=c("Down", "Not", "Up"))
gg1=ggbarplot(outTab, x="Pathway", y="t", fill = "Sig", color = "white",
              palette=c("green","black","#E41A1C"), sort.val = "asc", sort.by.groups = T,
              rotate=TRUE, legend="right", title=gene,
              xlab="Term", ylab="t value of GSVA score",  legend.title="Group", x.text.angle=60)
print(gg1)
dev.off()
