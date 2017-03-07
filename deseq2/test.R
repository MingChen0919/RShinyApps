library(DESeq2)
library("pasilla")
pasCts <- system.file("extdata", "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
# pasAnno <- system.file("extdata", "pasilla_sample_annotation.csv",
#                        package="pasilla", mustWork=TRUE)
countData <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
write.csv(countData, file='countData.csv', row.names = TRUE)
condition = c("treated", "treated", "treated", "untreated", "untreated", "untreated", "untreated")
type = c("single","paired","paired","single","single","paired","paired")
colData = data.frame(condition, type)
rownames(colData) = c("treated1","treated2","treated3",
                      "untreated1","untreated2","untreated3","untreated4")

## make sure rows in colData and columns in countData are in the same order
all(rownames(colData) %in% colnames(countData))
countData = countData[, rownames(colData)]
write.csv(colData, file='colData.csv', row.names = TRUE)


# construct a DESeqDataSet
dds = DESeqDataSetFromMatrix(countData = countData,
                             colData = colData,
                             design = ~ condition)
# pre-filtering (remove rows that have only 0 or 1 read)
dds
dds = dds[rowSums(counts(dds)) > 1, ]
dds

# set reference level
dds$condition = relevel(dds$condition, ref="untreated")

# differential expression analysis
dds = DESeq(dds)
res = results(dds)
res
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)

res05 = results(dds, alpha = 0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)


## MA-plot 
plotMA(res, main="DESeq2", ylim=c(-2,2))
df = data.frame(mean=res$baseMean,
                lfc = res$log2FoldChange,
                sig = ifelse(res$padj < 0.1, TRUE, FALSE))
plotMA(df, ylim=c(-2,2))


## plot counts
plotCounts(dds, gene=which.min(res$padj), intgroup = "condition")

countData['FBgn0039155', ]



##===== MA plot with plotly ========
library(plotly)
df = data.frame(gene_id = rownames(res),
                mean = res$baseMean,
                lfc = res$log2FoldChange,
                sig = ifelse(res$padj < 0.1, TRUE, FALSE),
                col = ifelse(res$padj < 0.1, "padj<0.1", "padj>=0.1"))
df$col = factor(df$col, ordered = FALSE)
df = subset(df, mean != 0)
df = subset(df, !is.na(sig))
# df = subset(df, sig == TRUE)
plot_ly(
  df, x = ~log(mean), y = ~lfc,
  # hover text:
  text = ~paste('Gene ID: ', gene_id, '<br>Normalized mean: ', mean, '<br>log fold change: ', lfc),
  color = ~col, colors = c("red", "blue")
  ) %>%
  layout(title = 'MA plot',
         xaxis = list(title="Log mean"),
         yaxis = list(title="log fold change"))

#option 2: ggplotly
p = ggplot(data = df, aes(x = log(mean), y = lfc, col = col)) +
  geom_point(aes(text = paste('Gene ID: ', gene_id, '<br>Normalized mean: ', mean, '<br>log fold change: ', lfc))) +
  xlab('Log base mean') + 
  ylab('Log fold change')
p = ggplotly(p)
p
##======== end of MA plot with plotly ===============



#========= plot count ===============
gene = 'FBgn0025111'
gene_counts = counts(dds[gene, ])
gene_counts
count = counts(dds[gene, ])[1,]
sample = colnames(gene_counts)
group = factor(colData[sample, 'condition'])
df_gene = data.frame(count, sample, group)
df_gene
plot_ly(
  type = 'scatter',
  df_gene, x = ~group, y = ~count, 
  color = ~group,
  text = ~paste('Sample: ', sample)
) %>% 
  layout(
    title = gene
  )


## option 2: ggplotly
p = ggplot(data=df_gene, aes(x=group, y=count, col=group)) +
  geom_jitter(aes(text=paste('Sample: ', sample)))
p = ggplotly(p)
p
#========= end of plot count ============




## heatmap of the count matrix
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("condition","type")]) 
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
# data
df
log2.norm.counts = as.data.frame(log2.norm.counts)
plot_ly(x=colnms, y=rownms, z = as.matrix(log2.norm.counts), 
        #key = as.matrix(log2.norm.counts), 
        type = "heatmap", source = "heatplot") %>%
  layout(xaxis = list(title = ""), 
         yaxis = list(title = ""))



library(DT)
dt_res = as.matrix(res)
datatable(dt_res, filter='top')

