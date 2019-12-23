setwd("C:/Users/Russell/OneDrive - UW/Linguistics/R/")

cts <- as.matrix(read.csv("bookCounts.tsv",sep="\t",row.names="Row",quote=""))
coldata <- read.csv("bookAnnotation.tsv",sep="\t", row.names=1)

cts <- as.matrix(read.csv("charCounts.tsv",sep="\t",row.names="Row",quote=""))
coldata <- read.csv("charAnnotation.tsv",sep="\t", row.names=1)

cts <- as.matrix(read.csv("bookCountsMystery.tsv",sep="\t",row.names="Row",quote=""))
coldata <- read.csv("bookAnnotationMystery.tsv",sep="\t", row.names=1)

cts <- as.matrix(read.csv("charCountsMystery.tsv",sep="\t",row.names="Row",quote=""))
coldata <- read.csv("charAnnotationMystery.tsv",sep="\t", row.names=1)


library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ gender + author + gender:author)


dds <- DESeq(dds, test="Wald", parallel=TRUE)




resultsNames(dds) # lists the coefficients
res <- results(dds, name="gender_M_vs_F")
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="gender_M_vs_F", type="apeglm")



write.csv(as.data.frame(res[order(res$pvalue),]), 
          file="gender_M_vs_F_char_Mystery_Interaction.csv")

write.csv(as.data.frame(res[order(res$pvalue),]), 
          file="gender_M_vs_F_char.csv")

write.csv(as.data.frame(res[order(res$pvalue),]), 
          file="gender_M_vs_F_char_Interaction_char.csv")

vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

write.csv(assay(vsd), 
          file="gender_M_vs_F_char_vsd.csv")









