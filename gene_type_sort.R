#################################################
#                                               #
# R script to sort genes according to gene type #
#            		 -Arindam Ghosh (19 March 2018) #
#                                               #
#################################################
# The input file is a *.tsv file as output of ballgown
# DEG analysis in R and gene type for each gene attached
# Contains columns: id, feature, fc, pval, qval, gene_name, Gene.type
#


# Load all required libraries 

library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)


# Load the input file, extract all gene types assigned

##### EDIT FILE NAME HERE #####
gene_list = read.table("gene_type_compare_gene_sigq_005.tsv", header= TRUE)
gene_type = as.data.frame(gene_list[,c("Gene.type")])
names(gene_type)[1] = "gene_type"

# Extract list of unique gene types
unique_gene_type = unique(gene_type, incomparables=FALSE)
rownames(unique_gene_type) = NULL


# New empty data frame to write results
result_write = data.frame()


# While loop to check the number of genes under ech gene type
# Write each set of gene types in separate file
# Write a file with the total counts under each gene type


i=1
while (i<= nrow(unique_gene_type) )
{

gene_type_list = subset(gene_list, Gene.type == unique_gene_type$gene_type[i])
##### EDIT FILE NAME HERE #####
write.table(gene_type_list, file=paste0(unique_gene_type$gene_type[i],"_gene_sigq_005.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
result_count = c(nrow(gene_type_list))
result_type = c(paste0(unique_gene_type$gene_type[i]))
result_write = rbind(result_write, data.frame(result_type, result_count))


i=i+1
}

result_write
##### EDIT FILE NAME HERE #####
write.table(result_write, "result_count.tsv", sep = "\t", quote = FALSE, row.names = FALSE, append =TRUE)
