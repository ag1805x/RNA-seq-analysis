#############################################################
#  Rscript for Differential Gene Expression Analysis using  #
#  DESeq2, EdgeR and Limma                                  #
#                        -Arindam Ghosh (16 November 2018)  #
#############################################################


library(edgeR)
library(DESeq2)
library(geneplotter)
library(gplots)
library(ggplot2)
library(pheatmap)
library(UpSetR)
library(RColorBrewer)


#Read count table and extract counts for Protein coding genes

#$$$ Change file name here
counts = read.table("PRJNA230824_allCounts_noOver_noMulti.tsv", row.names=c(1), header=T)
sprintf ("Total Genes")
dim(counts)
#$$$ Change path here
protein_coding = as.matrix(read.table("/home/bioinfo/Documents/Work/reference_genome/Human_84/Annotation/gff/protein_coding"))
sprintf ("Total Protein Coding Genes")
dim(protein_coding)
counts_pc <- counts[protein_coding,]
dim(counts_pc)



# Filtering out genes which donot have 1 cpm atleast in two samples
keep <- rowSums(cpm(counts_pc)>1) >1
counts_pc_filt <- counts_pc[keep, ]
sprintf ("Number of genes after filtering low count genes")
dim(counts_pc_filt)
#$$$ Change number of samples here
condition = c(rep("T",4), rep("C",2))
coldata <- data.frame(row.names=colnames(counts_pc_filt), condition)
sprintf ("Coldata")
coldata



# DESEQ2 Analysis
dds = DESeqDataSetFromMatrix(countData=counts_pc_filt, colData=coldata, design=~condition)
dds = DESeq(dds)
res = results(dds)
sprintf ("DESeq2 results")
dim(res)
#$$$ Change file name here (3)
write.table(res, "DS5_DESeq_result.tsv", sep="\t", col.names=NA)
res_clean = na.exclude(as.data.frame(res))
sprintf ("Number of genes after NA exclude (DESeq)")
dim(res_clean)
upreg = res_clean[(res_clean$log2FoldChange>1 & res_clean$padj<0.01),]
write.table(upreg, "DS5_DESeq_upreg.tsv", sep="\t", col.names=NA)
sprintf ("Number of upregulated genes (DESeq2)")
dim(upreg)
downreg = res_clean[(res_clean$log2FoldChange<(-1) & res_clean$padj<0.01),]
write.table(downreg, "DS5_DESeq_downreg.tsv", sep="\t", col.names=NA)
sprintf ("Number of downregulated genes (DESeq2)")
dim(downreg)



# Volcano Plot

#$$$ Change file name here 
png(filename="DS5_volcano.png")
with(res_clean, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano Plot", col="grey", xlim=c(-10,10)))
with(subset(res_clean,padj<0.01 & log2FoldChange>1),points(log2FoldChange, -log10(padj), pch=20, col="green"))
with(subset(res_clean,padj<0.01 & log2FoldChange<(-1)),points(log2FoldChange, -log10(padj), pch=20, col="red"))
abline(h=2, lty=2)
abline(v=-1, lty=2)
abline(v=1, lty=2)
dev.off()



# ECDF Plot

#$$$ Change file name here
png(filename="DS5_ecdf.png", height=1020, width=1020, units='px')
multiecdf(counts(dds, normalized=TRUE)[,], xlab="MeanCounts", xlim=c(0,1000))
dev.off()



# Plot Density

#$$$ Change file name here
png(filename="DS5_density.png", height=1020, width=1020, units='px')
multidensity(counts(dds, normalized=TRUE)[,], xlab="MeanCounts", xlim=c(0,1000))
dev.off()



# PLOT SAMPLE TO SAMPLE DISTANCE

#$$$ Change file name here (2)
rld = rlogTransformation(dds, blind=T)
png(filename="DS5_sample_heatmap1.png", height=1020, width=1020, units='px')
distRL = dist(t(assay(rld)))
mat=as.matrix(distRL)
hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13,13))
dev.off()
png(filename="DS5_sample_heatmap2.png", height=1020, width=1020, units='px')
pheatmap(mat, clustering_distance_rows=distRL, clustering_distance_cols=distRL, col=rev(hmcol))
dev.off()



# PCA

#$$$ Change file name here ()
png(filename="DS5_pca1.png")
pca = DESeq2::plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
p<-ggplot(pca,aes(x=PC1,y=PC2,color=group, label=row.names(pca) ))
p<-p+geom_point()
p
dev.off()

png(filename="DS5_pca2.png")
plotPCA(rld, intgroup=c("condition"))
dev.off()



# PLOT DISPERSION ESTIMATE

#$$$ Change file name here
png(filename="DS5_Dispersion.png", height=1020, width=1020, units='px')
plotDispEsts(dds)
dev.off()



# MA Plot

#$$$ Change file name here
png(filename="DS5_MA.png", height=1020, width=1020, units='px')
plotMA(res)
dev.off()



# Box-plot of Normalization

#$$$ Change file name here (3)

png(filename="DS5_boxplot_normalized.png",width = 1500, height= 1020, units="px")
#$$$ Change sample size
colors = c(rep("green",4),rep("blue",2))
boxplot(log2(counts(dds, normalized=TRUE)+1), col=colors, outline = FALSE, main="Box-plot of Normalized counts", xlab="Samples", ylab="log transformed normalized counts")
legend("topright", inset=0, title="Sample type", c("Pluripotent","Non-pluripotent"), fill=c("green","blue"), cex=0.8)
dev.off()

png(filename="DS5_boxplot_unnormalized.png",width = 1500, height= 1020, units="px")
boxplot(log2(counts_pc_filt+1), col=colors, outline = FALSE, main="Box-plot of Un-normalized counts", xlab="Samples", ylab="log transformed normalized counts")
legend("topright", inset=0, title="Sample type", c("Pluripotent","Non-pluripotent"), fill=c("green","blue"), cex=0.8)
dev.off()

png(filename="DS5_boxplot_normalized_rld.png",width = 1500, height= 1020, units="px")
#$$$ Change sample size
colors = c(rep("green",4),rep("blue",2))
boxplot(assay(rld), col=colors, outline = FALSE, main="Box-plot of Normalized counts", xlab="Samples", ylab="log transformed normalized counts")
legend("topright", inset=0, title="Sample type", c("Pluripotent","Non-pluripotent"), fill=c("green","blue"), cex=0.8)
dev.off()



#$$$ Change file name here
save (dds, file = 'DS5_dds.rda')



#   EDGER

#$$$ Change sample size
sample_info.edgeR <- factor(c(rep("T", 4), rep("C",2)))
sample_info.edgeR <- relevel(sample_info.edgeR, ref="C")
sprintf ("EdgeR sample info")
sample_info.edgeR

#$$$ Change file name here (4)
edgeR.DGElist <- DGEList(counts= counts_pc_filt, group=sample_info.edgeR)
sprintf ("Genes for edgeR")
dim(edgeR.DGElist$counts)

edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method="TMM")
edgeR.DGElist$samples
design <- model.matrix(~sample_info.edgeR)
sprintf ("edgeR design")
design
edgeR.DGElist <- estimateDisp(edgeR.DGElist, design)
edger_fit <- glmFit(edgeR.DGElist, design)
edger_lrt <- glmLRT(edger_fit)
DGE.results_edgeR <- topTags(edger_lrt, n=Inf, sort.by = "none", adjust.method = "BH")
write.table(DGE.results_edgeR$table, "DS5_EdgeR_result.tsv", sep="\t", col.names=NA)

edger.upreg = DGE.results_edgeR$table[(DGE.results_edgeR$table$logFC>1 & DGE.results_edgeR$table$FDR<0.01),]
write.table(edger.upreg, "DS5_EdgeR_upreg.tsv", sep="\t", col.names=NA)
sprintf ("Number of upregulated genes (edgeR)")
dim(edger.upreg)

edger.downreg = DGE.results_edgeR$table[(DGE.results_edgeR$table$logFC<(-1) & DGE.results_edgeR$table$FDR<0.01),]
write.table(edger.downreg, "DS5_EdgeR_downreg.tsv", sep="\t", col.names=NA)
sprintf ("Number of downregulated genes (edgeR)")
dim(edger.downreg)

#$$$ Change file name here
save (edger_lrt, file = 'DS5_edgeLRT.rda')



#  limma


rownames(design) <- colnames(edgeR.DGElist)
voomTransformed <- voom(edgeR.DGElist, design, plot=FALSE)
voom.fitted <- lmFit(voomTransformed, design = design)
voom.fitted <- eBayes(voom.fitted)
colnames(design)

#$$$ Change file name here (3)
DGE.results_limma <- topTable(voom.fitted, coef="sample_info.edgeRT", number=Inf, adjust.method="BH", sort.by="none")
write.table(DGE.results_limma, "DS5_limma_result.tsv", sep="\t", col.names=NA)

#$$$ Change file name here (2)
limma.upreg = DGE.results_limma[(DGE.results_limma$logFC>1 & DGE.results_limma$adj.P.Val<0.01),]
write.table(limma.upreg, "DS5_limma_upreg.tsv", sep="\t", col.names=NA)
sprintf ("Number of upregulated genes (limma)")
dim(limma.upreg)

limma.downreg = DGE.results_limma[(DGE.results_limma$logFC<(-1) & DGE.results_limma$adj.P.Val<0.01),]
write.table(limma.downreg, "DS5_limma_downreg.tsv", sep="\t", col.names=NA)
sprintf ("Number of downregulated genes (limma)")
dim(limma.downreg)

#$$$ Change file name here
save (voom.fitted, file = 'DS5_VoomFitted.rda')


# Results comparision

#$$$ Change file name here (2)
upreg_list <- list(DESeq2 = rownames(upreg), edgeR = rownames(edger.upreg), limma = rownames(limma.upreg))
upreg_gns <- UpSetR::fromList(upreg_list)
png(filename="DS5_upreg_comparision.png", height=720, width=720, units='px')
UpSetR::upset(upreg_gns, order.by="freq", text.scale = 2)
dev.off()

downreg_list <- list(DESeq2.downreg =rownames(downreg), limma.downreg=rownames(limma.downreg), edger.downreg=rownames(edger.downreg))
downreg_gns = fromList(downreg_list)
png(filename="DS5_downreg_comparision.png", height=720, width=720, units='px')
upset(downreg_gns, order.by="freq", text.scale = 2)
dev.off()

