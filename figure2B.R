install.packages("pheatmap")
install.packages("reshape2")
install.packages("ggpubr")

rm(list=ls())

library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

expFile="normalize_GSE23649.txt"


#
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

data=data[rowMeans(data)>0,]
exp=data

#
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))
conNum=length(Type[Type=="Control"])
treatNum=length(Type[Type=="Treat"])

#
sigVec=c()
sigGeneVec=c()
diffTab=data.frame()
for(i in row.names(data)){
  test=wilcox.test(data[i,] ~ Type)
  pvalue=test$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  conMean=mean(data[i,1:conNum])
  treatMean=mean(data[i,(conNum+1):ncol(data)])
  logFC=treatMean-conMean
  if(pvalue<0.05){
    sigVec=c(sigVec, paste0(i, Sig))
    sigGeneVec=c(sigGeneVec, i)
    if(logFC>0){diffTab=rbind(diffTab, cbind(Gene=i,conMean,treatMean,pvalue,logFC,Type="Up"))}
    if(logFC<0){diffTab=rbind(diffTab, cbind(Gene=i,conMean,treatMean,pvalue,logFC,Type="Down"))}
  }
}

pvalue=diffTab[,"pvalue"]
fdr=p.adjust(as.numeric(as.vector(pvalue)),method="fdr")
diffTab=cbind(diffTab,fdr=fdr)

#
write.table(diffTab, file="diff.txt", sep="\t", quote=F, row.names=F)




###
fdrFilter=0.05                                                    #fdr?ٽ?ֵ
logFCfilter=1                                                    #logFC?ٽ?ֵ  ???޸İ??ղ????????????޸?

diffTab=diffTab[( abs(as.numeric(as.vector(diffTab$logFC)))>logFCfilter & as.numeric(as.vector(diffTab$fdr))<fdrFilter),]
write.table(diffTab,file="diff_significant.xls",sep="\t",row.names=F,quote=F)


#a=data[sigGeneVec,]
outTab=rbind(ID=colnames(data), data)
write.table(outTab, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)
row.names(data)=sigVec

#?Բ?es(Type)=colnames(data)
Type=as.data.frame(Type)

Type <- cbind(rownames(Type),Type)  

Type=Type[order(Type$Type),]
Type$AA<- "A"
Type1<-  as.data.frame(Type[,2])  
class(Type)
rownames(Type1) <- Type[,1] 
colnames(Type1)[1]<- "Type"
Type<-Type1

data<-t(data)
data <- cbind(rownames(data),data) 
Type <- cbind(rownames(Type),Type) 
colnames(Type)

colnames(data)[1]<-"rownames(Type)"
data<- merge(Type,data, by="rownames(Type)")
data=data[order(data$Type),]
rownames(data) <- data[,1] 
data<-data[,-c(1,2)]
data<-t(data)
class(data)
data<-as.data.frame(data)
Type<-Type1
write.table(Type,"Type.txt",sep = "\t",quote = F)

write.table(data,"data.txt",sep = "\t",quote = F)
Type= read.table("Type.txt",header=T,sep="\t",row.names=1,check.names=F)

data= read.table("data.txt",header=T,sep="\t",row.names=1,check.names=F)



pdf(file="heatmap.pdf", width=9, height=5)
pheatmap(data,
         annotation =Type,
         color = colorRampPalette(c(rep("green",3), "white", rep("red",3)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=F,
         fontsize=6,
         fontsize_row=7,
         fontsize_col=6)
dev.off()
