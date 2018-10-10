#Load the required libraries
library(DESeq2)
library(geneplotter)
library(dplyr)
library(genefilter)
library(RColorBrewer)
library(stringr)
library(gplots)
library(LSD)
library(ggplot2)
library(pheatmap)


#Read in the count data
count = read.table("count.tsv", header=T, row.names=c(1))
sprintf ("Total number of genes read in:")
dim(count)

#Filter genes that have zero counts across all samples
count_filt = count[rowSums(count)>0,]
sprintf ("Genes ofter low count filter:")
dim(count_filt)

#Define conditions for the samples
condition = c(rep("c",3), rep("t",3))
coldata = data.frame(colnames(count_filt), condition)
coldata

#Create DESeq object for analysis with DESeq2
dds = DESeqDataSetFromMatrix(countData=count_filt, colData=coldata, design= ~ condition)

#DESeq2 analysis
dds = DESeq(dds)

#Generate results
res = results(dds)
write.table(res, "DESeq_result.tsv", sep="\t", col.names=NA)
res_clean = na.exclude(as.data.frame(res))
sprintf ("Genes after NA EXCLUDE")
dim(res_clean)


upreg = res_clean[(res_clean$log2FoldChange>1 & res_clean$padj<0.01),]
write.table(upreg, "upreg.tsv", sep="\t", col.names=NA)
sprintf ("UPREG Genes")
dim(upreg)


downreg = res_clean[(res_clean$log2FoldChange<(-1) & res_clean$padj<0.01),]
write.table(downreg, "downreg.tsv", sep="\t", col.names=NA)
sprintf ("DOWNREG Genes")
dim(downreg)


#VOLCANO PLOT
pdf("volcano.pdf")
with(res_clean, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano Plot", col="grey", xlim=c(-10,10)))
with(subset(res_clean,padj<0.01 & log2FoldChange>1),points(log2FoldChange, -log10(padj), pch=20, col="green"))
with(subset(res_clean,padj<0.01 & log2FoldChange<(-1)),points(log2FoldChange, -log10(padj), pch=20, col="red"))
#with(top,points(log2FoldChange, -log10(padj), pch=20, col="blue"))
abline(h=2, lty=2)
abline(v=-1, lty=2)
abline(v=1, lty=2)
dev.off()



# PLOT ECDF
pdf("ecdf.pdf")
multiecdf(counts(dds, normalized=TRUE)[,], xlab="MeanCounts", xlim=c(0,1000))
dev.off()

#PLOT DENSITY
pdf("density.pdf")
multidensity(counts(dds, normalized=TRUE)[,], xlab="MeanCounts", xlim=c(0,1000))
dev.off()

# PLOT SAMPLE TO SAMPLE DISTANCE
rld = rlogTransformation(dds, blind=T)
pdf("sample_heatmap1.pdf")
distRL = dist(t(assay(rld)))
mat=as.matrix(distRL)
hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13,13))
dev.off()
pdf("sample_heatmap2.pdf")
pheatmap(mat, clustering_distance_rows=distRL, clustering_distance_cols=distRL, col=hmcol)
dev.off()

# PLOT PCA
pdf("pca.pdf")
pca = DESeq2::plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
p<-ggplot(pca,aes(x=PC1,y=PC2,color=group, label=row.names(pca) ))
p<-p+geom_point()
p
dev.off()

# PLOT DISPERSION ESTIMATE
pdf("Dispersion.pdf")
plotDispEsts(dds)
dev.off()


# MA Plot
pdf("MA.pdf")
plotMA(res)
dev.off()

# Box-plot of Normalization
pdf("boxplot_normalized.pdf")

colors = c(rep("blue",3), rep("green",3))
boxplot(log2(counts(dds, normalized=TRUE)+1), col=colors, outline = FALSE, main="Box-plot of Normalized counts", xlab="Samples", ylab="log transformed normalized counts")
legend("topright", inset=0, title="Sample type", c("Pluripotent","Non-pluripotent"), fill=c("green","blue"), cex=0.8)
dev.off()
pdf("boxplot_unnormalized.pdf")
boxplot(log2(count+1), col=colors, outline = FALSE, main="Box-plot of Un-normalized counts", xlab="Samples", ylab="log transformed normalized counts")
legend("topright", inset=0, title="Sample type", c("Pluripotent","Non-pluripotent"), fill=c("green","blue"), cex=0.8)
dev.off()


save (dds, file = 'dds.rda')





