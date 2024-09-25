install.packages("rms")
install.packages("rmda")


#
library(rms)
library(rmda)

inputFile="normalize_GSE23649.txt"          
geneFile="interGenes.txt"    


#
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
row.names(data)=gsub("-", "_", row.names(data))

#eRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

#??È¡a=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type=group)
paste(colnames(data), collapse="+")

#???Ýst=datadist(rt)
options(datadist="ddist")

#????Model=lrm(Type~ SLC7A5+SLC1A5+SLC43A2, data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
	fun.at=c(0.0001,0.1,0.3,0.5,0.7,0.9,0.99),
	lp=F, funlabel="Risk of Disease")
#????("Nomo.pdf", width=8, height=6)
plot(nomo)
dev.off()


#????i=calibrate(lrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=5.5, height=5.5)
plot(cali,
	xlab="Predicted probability",
	ylab="Actual probability", sub=F)
dev.off()


#???ÆType=ifelse(rt$Type=="Control", 0, 1)
dc=decision_curve(Type ~ SLC7A5+SLC1A5+SLC43A2, data=rt, 
	family = binomial(link ='logit'),
	thresholds= seq(0,1,by = 0.01),
	confidence.intervals = 0.95)
#????(file="DCA.pdf", width=5.5, height=5.5)
plot_decision_curve(dc,
	curve.names="Model",
	xlab="Threshold probability",
	cost.benefit.axis=T,
	col="green",
	confidence.intervals=FALSE,
	standardize=FALSE)
dev.off()



#????(file="clinical_impact.pdf", width=6, height=6)
plot_clinical_impact(dc,
                     confidence.intervals=T,
                     col = c("green", "red"))
dev.off()

####