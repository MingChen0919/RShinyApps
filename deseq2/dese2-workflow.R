## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("tximport")
ainstall.packages('readr')
library(tximport)
library(readr)
library(DESeq2)

## import count matrix data
count_matrix = read.csv()

## pre-processing data
#  remove rows that have only 0 or 1 read.
dds = DESeqDataSetFromMatrix()
dds = dds[rowSums(ounts(dds)) > 1, ]
# factor levels
dds$condition = relevel(dds$condition, ref="untreated")
# collapsing technical replicates (NOT biological replicates)
collapseReplicates()

##====columns of the count matrix and rows of the column data have to be in the same order!


## differential expression analysis
DESeq()