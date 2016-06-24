
#Install the DESeq, gplots, and pasilla libraries if needed, then load.


#system.file() used to get the path to the directory containing the count table.

datafile <- system.file( "extdata/pasilla_gene_counts.tsv", package="pasilla" )

#Load file using read.table: creates a data frame from it. Store as count.table. Header = TRUE if 
#the first row of your file contains the names of variables.

count.table <- read.table(datafile, header=TRUE, row.names=1)

#See the first few rows of data.
head(count.table)

#Optional: delete any genes that were never detected (equal to zero in all conditions).
#count.table <- count.table[rowSums(count.table) > 0,]

#Metadata. Columns will be condition or libType. Rows correspond to the 7 (untreated or treated) samples.
pasillaDesign = data.frame(
  row.names = colnames( count.table ), 
  condition = c( "untreated", "untreated", "untreated",
                 "untreated", "treated", "treated", "treated" ), 
  libType = c( "single-end", "single-end", "paired-end",
               "paired-end", "single-end", "paired-end", "paired-end" ) )

#This dataset contains RNA-Seq counts for RNAi treated vs untreated cells, as well as sequencing data
#using both single-end sequencing (reading fragments from only one end to the other) and paired-end 
#sequencing (reading form both ends).

pairedSamples = pasillaDesign$libType == "paired-end"
countTable = count.table[ , pairedSamples ]
condition = pasillaDesign$condition[ pairedSamples ]

head(countTable)

cds = newCountDataSet( countTable, condition )

#Normalize the expression values of each treatment by dividing each column with its own size factor using
#the estimateSizeFactors function. This estimates the size factor by first taking each column and dividing
#by the geometric mean of the rows. The median of these ratios is used as the size factor for this column.

cds = estimateSizeFactors( cds )
sizeFactors( cds )
#This divides each column by the size factor, which makes them comparable.
head( counts( cds, normalized=TRUE ) )


#Estimate dispersion parameter
cds = estimateDispersions( cds )

#Plotting dispersion estimates. 
plotDispEsts( cds )

#See whether there is differential expression between untreated and treated. 
#Output is a data.frame. 
res = nbinomTest( cds, "untreated", "treated" )
head(res)

#To order by p-vals (decreasing)
res <- res[order(res$padj),]
head(res)

#Plot log2fold changes against mean normalised counts for untreated vs treated.
plotMA(res, col = ifelse(res$padj >=0.01, "black", "violet"))
abline(h=c(-1:1), col="red")

#Select gene names based on FDR (1%)
gene.kept <- res$id[res$padj <= 0.01 & !is.na(res$padj)]

#Create a count data set with multiple factors
cdsFull = newCountDataSet(count.table, pasillaDesign)

#Estimating size factors
cdsFull = estimateSizeFactors( cdsFull )
#Estimating dispersions
cdsFull = estimateDispersions( cdsFull )

#Variance stabilizing transformation
cdsBlind = estimateDispersions( cds, method="blind" )
vsd = varianceStabilizingTransformation( cdsBlind )

#Heatmap of the count table from the variance stabilisation transformed data for the 30 most highly expressed genes.
cdsFullBlind = estimateDispersions( cdsFull, method = "blind" )
vsdFull = varianceStabilizingTransformation( cdsFullBlind )

select = order(rowMeans(counts(cdsFull)), decreasing=TRUE)[1:30]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(10, 6))

#Sample clustering
distances = dist( t( exprs(vsdFull) ) )

mat = as.matrix( distances )
rownames(mat) = colnames(mat) = with(pData(cdsFullBlind), paste(condition, libType, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

#PCA plot of the samples
print(plotPCA(vsdFull, intgroup=c("condition", "libType")))

