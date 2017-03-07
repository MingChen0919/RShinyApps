library(shiny)
library(DT)
library(plotly)


function(input, output) {
  # countData
  countData = reactive({
    if(is.null(input$countDataInput)){
      return(NULL)
    }
    countData = read.csv(input$countDataInput$datapath, header=TRUE, row.names=1)
  })
  
  # colData
  colData = reactive({
    if(is.null(input$colDataInput)){
      return(NULL)
    }
    
    read.csv(input$colDataInput$datapath, header=TRUE, row.names=1)
  })
  
  # construct a DESeqDataSet
  dds = reactive({
    dds = DESeqDataSetFromMatrix(countData = countData(),
                                 colData = colData(),
                                 design = ~ condition)
    # pre-filtering (remove rows that have only 0 or 1 read)
    dds = dds[rowSums(counts(dds)) > 1, ]
    # set reference level
    dds$condition = relevel(dds$condition, ref="untreated")
    # differential expression analysis
    DESeq(dds)
  })


  #>>>render MA plot
  output$MA_plot = renderPlotly({
    res = results(dds())
    # df for plot
    df = data.frame(gene_id = rownames(res),
                    mean = res$baseMean,
                    lfc = res$log2FoldChange,
                    sig = ifelse(res$padj < 0.1, TRUE, FALSE),
                    col = ifelse(res$padj < 0.1, "padj<0.1", "padj>=0.1"),
                    stringsAsFactors = FALSE)
    p = ggplot(data = df, aes(x = log(mean), y = lfc, col = col, key = gene_id)) +
      geom_point(aes(text = paste('Gene ID: ', gene_id, '<br>Normalized mean: ', mean, '<br>log fold change: ', lfc))) +
      xlab('Log base mean') + 
      ylab('Log fold change')
    ggplotly(p) %>% layout(dragmode = 'select',
                           title = "MA plot")
    #=====!!! plot_ly() reports wrong keys======
    # plot_ly(
    #   df, x = ~log(mean), y = ~lfc, key = ~gene_id,
    #   # hover text:
    #   text = ~paste('Gene ID: ', gene_id, '<br>Base mean: ', mean, '<br>log fold change: ', lfc),
    #   color = ~col, colors = c("red", "blue")
    # ) %>%
    #   layout(title = 'MA plot',
    #          xaxis = list(title="Log mean"),
    #          yaxis = list(title="Log fold change"),
    #          dragmode = "select")
  })
  
  #>>>render countData table
  output$countDataTable = renderDataTable({
    datatable(countData(), filter='top', caption='Table 1: gene expression count matrix')
  })
  
  #>>>render colData table
  output$colDataTable = renderDataTable({
    datatable(colData(), filter='top', caption='Table 2: gene expression treatment information')
  })
  
  #>>>render click and drag data
  output$brush = renderDataTable({
    d = event_data("plotly_selected")
    res = results(dds())[d$key, ]
    res = as.data.frame(res)
    # if (is.null(d)) "Click and drag events appear here" else datatable(res)
    datatable(res)
  })
  
  #>>>render click
  output$click = renderPlotly({
    d = event_data("plotly_click")
    dds = dds()
    dds_count = counts(dds[d$key, ])
    count = dds_count[d$key,]
    sample = colnames(dds_count)
    group = colData()[sample, 'condition']
    df = data.frame(count,sample, group)
    # ggplotly
    p = ggplot(data=df, aes(x=group, y=count, col=group)) +
      geom_jitter(aes(text=paste('Sample: ', sample)))
    ggplotly(p) %>% layout(dragmode = 'select',
                           title = d$key,
                           xaxis = list(title="Group"),
                           yaxis = list(title="Count"))

  })
}