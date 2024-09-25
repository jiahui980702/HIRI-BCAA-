if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
 
install.packages("pheatmap")


library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)




#Read the presentation data file
rt=read.table("CRGexp.txt", header=T, sep="\t", check.names=F)

rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
exp=data

#Sample grouping information was extracted
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))

#Difference gene analysis
sigVec=c()
sigGeneVec=c()
for(i in row.names(data)){
  test=wilcox.test(data[i,] ~ Type)
  pvalue=test$p.value
  Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
  if(pvalue<0.05){
    sigVec=c(sigVec, paste0(i, Sig))
    sigGeneVec=c(sigGeneVec, i)}
}
#Output the expression of differential genesa=data[sigGeneVec,]
outTab=rbind(ID=colnames(data), data)
write.table(outTab, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)
row.names(data)=sigVec

#?Բ?The differential genes were visualized and heat mappedes(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=7, height=4.5)
pheatmap(data,
         annotation=Type,
         color = colorRampPalette(c(rep("green",3), "white", rep("red",3)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=7,
         fontsize_row=7,
         fontsize_col=7)
dev.off()

#?ѱ?Convert the presentation data to a ggplot2 input file=as.data.frame(t(exp))
exp=cbind(exp, Type=Type)
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#????Box plottinggboxplot(data, x="Gene", y="Expression", color = "Type", 
            xlab="",
            ylab="Gene expression",
            legend.title="Type",
            palette = c("green", "red"),
            add="point",
            width=0.8)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#????Output box diagram(file="boxplot.pdf", width=8, height=6)
print(p1)
dev.off()


##
