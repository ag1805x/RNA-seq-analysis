##############################################
#                                            #
#  Rscript to combine output of htseq-count  #
#                       -Arindam Ghosh       #
#                        (7 June 2018)       #
#                                            #
##############################################




#Get list of all files in directory

files = list.files(path=".", pattern="SRR")

#Merge the first to files and store
file1 = read.table(files[1], col.names=c("gene_id",files[1]))
file2 = read.table(files[2], col.names=c("gene_id",files[2]))
out.file = merge (file1, file2, by=c("gene_id"))

#For loop to merge contents of remaining files

for(i in 3:length(files))
{
file = read.table(files[i],col.names=c("gene_id",files[i]))
out.file <- merge(out.file, file, by=c("gene_id"))
}

write.table(out.file, file = "htseq_all_sample_count.tsv",sep="\t", row.names = FALSE)
