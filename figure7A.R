install.packages("pheatmap")


#
library(pheatmap)

# 

#
inputFile <- "ssGSEA.result.txt"  # ssGSEA
rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

#
con <- grepl("_Control", colnames(rt), ignore.case = TRUE)
treat <- grepl("_treat", colnames(rt), ignore.case = TRUE)
conData <- rt[, con]
treatData <- rt[, treat]
conNum <- ncol(conData)
treatNum <- ncol(treatData)
data <- cbind(conData, treatData)

#
data_matrix <- as.matrix(data)

#
Type <- c(rep("Control", conNum), rep("Treat", treatNum))
names(Type) <- colnames(data_matrix)
Type <- as.data.frame(Type)

#
pdf(file = "heatmap.pdf", width = 8, height = 5)
pheatmap(data_matrix, 
         annotation_col = Type,  # ?Զ???????ɫ
         color = colorRampPalette(c(rep("green", 2), "white", rep("red", 2)))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         scale = "row",
         show_colnames = FALSE,
         show_rownames = TRUE,
         fontsize = 7,
         fontsize_row = 7,
         fontsize_col = 7)
dev.off()
