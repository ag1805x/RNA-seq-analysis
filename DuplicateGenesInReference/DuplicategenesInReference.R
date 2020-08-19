# Identifying genes with same name in the Ensembl Reference Genome Annotation
# GeneAnnotation.csv based on GRCh38.p5

annot <- read.csv("GeneAnnotation.csv", stringsAsFactors=F, header=T)
temp <- as.data.frame(table(annot$AssociatedGeneName), stringsAsFactors=F)
colnames(temp)[1] <- "AssociatedGeneName"
annot.duplicate <- annot[annot$AssociatedGeneName %in% temp[temp$Freq>1,]$AssociatedGeneName, ]
annot.unique <- annot[annot$AssociatedGeneName %in% temp[temp$Freq==1,]$AssociatedGeneName, ]
annot.others <- annot[annot$AssociatedGeneName %in% temp[temp$Freq<1,]$AssociatedGeneName, ]
write.table(annot.duplicate, "DuplicateGenes.tsv", sep="\t", quote=F, row.names=F)
