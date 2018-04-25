#################################################
#                                               #
# R script to create adjacency matrix from gene #
# expression values to be used for creating     #
# network                                       #
#      -Arindam Ghosh (25 April 2018)           #
#                                               #
#################################################



#Read expression data from csv file
#Assumes file contains gene names/id in column 1 followed by expression values in next columns
exp_data = read.csv("data.csv", header = TRUE, row.names= c(1))


#Calculate Pearsons correlation between all genes
corr = cor (t(exp_data), method="pearson")


#Keep backup of correlation values
adj_matrx = corr

#Create adjacency matrix based on threshold 0.9

for(row in 1:nrow(adj_matrx)) 
{
for(col in 1:ncol(adj_matrx)) 
{
if (adj_matrx[row, col] > 0.9 ||   adj_matrx[row, col] < -0.9)    
{
adj_matrx[row, col] = 1
}
else
{
adj_matrx[row, col] = 0
}
}
}


#Write adjacency matrix to file
write.csv(adj_matrx, "adj_matrx.csv", row.names = TRUE)

