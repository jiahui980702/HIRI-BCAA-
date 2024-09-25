install.packages("ggplot2")
install.packages("ggrepel")


library(dplyr)
library(ggplot2)
library(ggrepel)


logFCfilter=0.5              
pvalue.Filter=0.05      
inputFile="diff.txt"    


#
rt = read.table(inputFile, header=T, sep="\t", check.names=F)
colnames(rt)
#
Sig=ifelse((rt$pvalue<pvalue.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")

#
rt = mutate(rt, Sig=Sig)
p = ggplot(rt, aes(logFC, -log10(pvalue)))+
  geom_point(aes(col=Sig))+
  scale_color_manual(values=c("#377EB8", "black","#E41A1C"))+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
#
p1=p+geom_label_repel(data=filter(rt, ((rt$pvalue<pvalue.Filter) & (abs(rt$logFC)>logFCfilter))),
                      box.padding=0.1, point.padding=0.1, min.segment.length=0.05,
                      size=1.8, aes(label=Gene)) + theme_bw()
#
pdf(file="vol.pdf", width=7, height=6.1)
print(p1)
dev.off()
