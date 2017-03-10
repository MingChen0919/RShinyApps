library(shiny)
library(DT)
library(plotly)


function(input, output, session) {
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
  
  ##----- display upload data in tables -------------
  #>>>render countData table
  output$countDataTable = renderDataTable({
    datatable(countData(), filter='top', caption='Table 1: gene expression count matrix')
  })
  
  #>>>render colData table
  output$colDataTable = renderDataTable({
    datatable(colData(), filter='top', caption='Table 2: gene expression treatment information')
  })
  
  
  #===== Wald Test ===============
  # add features to design formula
  output$waldDesign = renderUI({
    features = colnames(colData())
    names(features) = features
    selectizeInput('waldDesign', label='Add factors to model',
                   multiple=TRUE,
                   choices=features)
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
  # wald test
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
  
  
  #------ DESeq analysis complete check --------
  output$deseq_complete_check = reactive({
    if(is.na(ddsWald())) return(FALSE) else return(TRUE)
  })
  outputOptions(output, 'deseq_complete_check', suspendWhenHidden=FALSE)
  #---------------------------------------------
  
  #======== user input UI results exploration ======
  #>>>render factor selection UI
  output$waldContrastFactor = renderUI({
    factor = input$waldDesign
    selectInput('waldContrastFactor', label='Contrast Factor',
                choices=input$waldDesign)
  })
  #>>>render control level from selected factor
  output$waldControlLevel = renderUI({
    selectedFactor = input$waldContrastFactor
    allLevels = colData()[, selectedFactor]
    selectInput('waldControlLevel', label='Contrast Level',
                choices=allLevels)
  })
  
  #---control level selection complete check---
  output$waldControlLevelComplete = reactive({
    if(is.na(input$waldControlLevel)) return(FALSE) else return(TRUE)
  })
  outputOptions(output, 'waldControlLevelComplete', suspendWhenHidden=FALSE)
  
  #>>>render treated level from selected factor
  output$waldTreatedLevel = renderUI({
    selectedFactor = input$waldContrastFactor
    allLevels = colData()[, selectedFactor]
    # remove element that has been selected as control level
    allLevels = allLevels[allLevels != input$waldControlLevel]
    selectInput('waldTreatedLevel', label='Treated Level',
                choices=allLevels)
  })
  ##================================================
  
  
  ##======= render results after user input ========
  ## results(dds, contrast=c("condition","B","A"))
  contrastTable = reactive({
    contrastFactor = input$waldContrastFactor
    controlLevel = input$waldControlLevel
    treatedLevel = input$waldTreatedLevel
    res = results(ddsWald(), contrast=c(contrastFactor, controlLevel, treatedLevel))
    # res = as.matrix(res)
  })
  
  output$waldResultTable = renderDataTable({
    caption = paste0('Gene expression comparison: ', 
                     input$waldControlLevel, ' vs ', input$waldTreatedLevel)
    df = as.matrix(contrastTable())
    datatable(df, filter="top", caption=caption)
  })
  ##================================================
  
  
  
  ##========== Resualt Visualization ===============
  ##
  ##-----------MA plot -----------------------------
  dfMAplot = reactive({
    res = contrastTable()
    # df for plot
    df = data.frame(gene_id = rownames(res),
                    mean = res$baseMean,
                    lfc = res$log2FoldChange,
                    padj = res$padj,
                    # sig = ifelse(res$padj < 0.1, TRUE, FALSE),
                    col = ifelse(res$padj < 0.1, "padj<0.1", "padj>=0.1"),
                    stringsAsFactors = FALSE)
    rownames(df) = df$gene_id
    return(df)
  })
  # output$dfMAplot = renderPrint({
  #   dfMAplot()
  # })
  #>>>render MA plot
  output$MA_plot = renderPlotly({
    p = ggplot(data = dfMAplot(), aes(x = log(mean), y = lfc, col = col, key = gene_id)) +
      geom_point(aes(text = paste('Gene ID: ', gene_id, '<br>Normalized mean: ', mean, '<br>log fold change: ', lfc))) +
      xlab('Log base mean') +
      ylab('Log fold change')
    ggplotly(p) %>% layout(dragmode = 'select',
                           title = "MA plot")
  })
  
  # #>>>render click and drag data
  output$brush = renderDataTable({
    d = event_data("plotly_selected")
    if(is.na(d$key)) {
      return(NULL)
    } else {
      # res = subset(dfMAplot(), gene_id == d$key)
      df = dfMAplot()
      # res = df[df$gene_id == d$key, ]
      res = df[d$key, ]
      datatable(res)
    }
  })

  #>>>render click
  output$click = renderPlotly({
    d = event_data("plotly_click")
    if(is.null(d$key)) return(NULL)
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
  
  
  ##-------- heatmap --------------
  ## heatmap data
  output$heatmap = renderPlotly({
    cnt = counts(ddsWald())
    x = colnames(cnt)
    y = rownames(cnt)
    z = as.matrix(cnt)
    
    plot_ly(x=x, y=y, z=z, key=y,
            type="heatmap", source="heatplot") %>%
      layout(xaxis = list(title = ""), 
             yaxis = list(title = ""))
  })
  ## gene correlation from heatmap click
  # output$heatmapClick = renderPlotly({
  #   s <- event_data("plotly_click", source = "heatplot")
  #   if (length(s)) {
  #     vars <- c(s[["x"]], s[["y"]])
  #     d <- setNames(mtcars[vars], c("x", "y"))
  #     yhat <- fitted(lm(y ~ x, data = d))
  #     plot_ly(d, x = ~x) %>%
  #       add_markers(y = ~y) %>%
  #       add_lines(y = ~yhat) %>%
  #       layout(xaxis = list(title = s[["x"]]), 
  #              yaxis = list(title = s[["y"]]), 
  #              showlegend = FALSE)
  #   } else {
  #     plotly_empty()
  #   }
  # })
  

}