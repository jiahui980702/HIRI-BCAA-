install.packages("VennDiagram")


library(VennDiagram)      #
lassoFile="LASSO.gene.txt"      #lasso
RFFile="RFGenes.txt"      #
svmFile="SVM-RFE.gene.txt"      #


geneList=list()

#
rt=read.table(lassoFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              
geneNames=gsub("^ | $","",geneNames)     
uniqGene=unique(geneNames)               
geneList[["LASSO"]]=uniqGene             

#
rt=read.table(RFFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              
geneNames=gsub("^ | $","",geneNames)     
uniqGene=unique(geneNames)              
geneList[["RF"]]=uniqGene              


# 
# 
rt=read.table(svmFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])            
geneNames=gsub("^ | $","",geneNames)    
uniqGene=unique(geneNames)            
geneList[["SVM"]]=uniqGene            



#
venn.plot=venn.diagram(geneList,filename=NULL,fill=c("#984ea3", "red", "green"),
                       scaled=FALSE,cat.col = c("#984ea3", "red", "green"),cat.cex = 1.2)
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#
interGenes=Reduce(intersect,geneList)
write.table(file="interGenes.txt", interGenes, sep="\t", quote=F, col.names=F, row.names=F)

