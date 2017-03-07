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
  
  # file uploaded confirmation
  output$file_upload_check = reactive({
    if(is.null(input$countDataInput) | is.null(input$colDataInput)) return(NULL) else return(TRUE)
  })
  outputOptions(output, 'file_upload_check', suspendWhenHidden=FALSE)
  
  
  #===== Wald Test ===============
  # add features to design formula
  output$waldDesign = renderUI({
    features = colnames(colData())
    names(features) = features
    selectizeInput('waldDesign', label='Add factors to model',
                   multiple=TRUE,
                   choices=features)
  })
  # select contrast factors
  output$waldContrastFactor = renderUI({
    factor = input$waldDesign
    selectInput('waldContrastFactor', label='Contrast Factor',
                choices=input$waldDesign)
  })
  
  #===== LRT Test ===============
  # features in full model
  output$lrtFullModel = renderUI({
    features = colnames(colData())
    names(features) = features
    selectizeInput('lrtFullModel', label='Factors included in full model',
                   multiple=TRUE,
                   choices=features)
  })
  
  output$lrtReducedModel = renderUI({
    features = colnames(colData())
    names(features) = features
    selectizeInput('lrtReducedModel', label='Factors included in reduced model',
                   multiple=TRUE,
                   choices=features)
  })
  
  
  #========= Run DESeq Analysis ===========
  ddsWald = eventReactive(input$runWald, {
    countData = countData()
    colData = colData()
    design = as.formula(paste0('~', paste(input$waldDesign, collapse='+')))
    dds = DESeqDataSetFromMatrix(countData = countData,
                                 colData = colData,
                                 design = design)
    # pre-filtering (remove rows that have only 0 or 1 read)
    dds = dds[rowSums(counts(dds)) > 1, ]
    # set reference level
    # dds$condition = relevel(dds$condition, ref="untreated")
    # differential expression analysis
    DESeq(dds)
  })
  
  output$ddsWaldPrint = renderPrint({
    countData = countData()
    colData = colData()
    design = as.formula(paste0('~', paste(input$waldDesign, collapse='+')))
    dds = DESeqDataSetFromMatrix(countData = countData,
                                 colData = colData,
                                 design = design)
    # pre-filtering (remove rows that have only 0 or 1 read)
    dds = dds[rowSums(counts(dds)) > 1, ]
    # set reference level
    # dds$condition = relevel(dds$condition, ref="untreated")
    # differential expression analysis
    dds = DESeq(dds)
    res = results(dds)
    res
  })
  
  # render results as a datatable
  output$resTableWald = renderDataTable({
    res = results(ddsWald())
    datatable(as.matrix(res))
  })

  
  #>>>render countData table
  output$countDataTable = renderDataTable({
    datatable(countData(), filter='top', caption='Table 1: gene expression count matrix')
  })
  
  #>>>render colData table
  output$colDataTable = renderDataTable({
    datatable(colData(), filter='top', caption='Table 2: gene expression treatment information')
  })
  
  #>>>render MA plot
  output$MA_plot = renderPlotly({
    res = results(ddsWald())
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
  })
  
  #>>>render click and drag data
  output$brush = renderDataTable({
    d = event_data("plotly_selected")
    res = results(ddsWald())[d$key, ]
    res = as.data.frame(res)
    # if (is.null(d)) "Click and drag events appear here" else datatable(res)
    datatable(res)
  })
  
  #>>>render click
  output$click = renderPlotly({
    d = event_data("plotly_click")
    if(is.null(d)) return(NULL)
    dds = ddsWald()
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