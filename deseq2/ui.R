library(shiny)
library(plotly)
library(DESeq2)
library(DT)



fluidPage(
  sidebarLayout(
    sidebarPanel(
      #---data input---
      tags$h3('Data Input'),
      tags$hr(),
      fileInput('countDataInput', label='Upload count data'),
      fileInput('colDataInput', label='Upload column data'),
      #---DESeq analysis conditional on complete data upload
      conditionalPanel("output.file_upload_check",
                       tags$h3('DESeq Analysis'),
                       tags$hr(),
                       selectInput('testMethod', label='Wald or LRT test', 
                                   choices=c('Wald' = 'wald',
                                             'LRT' = 'lrt'),
                                   selected='wald'),
                       conditionalPanel("input.testMethod == 'wald'",
                                        uiOutput('waldDesign'),
                                        actionButton('runWald', label = 'Run')
                                        ),
                       conditionalPanel("input.testMethod == 'lrt'",
                                        uiOutput('lrtFullModel'),
                                        uiOutput('lrtReducedModel'),
                                        actionButton('runLRT', label = 'Run')
                                        )
                       ),
      #---DESeq results exploration conditional on complete DESeq analysis
      conditionalPanel('output.deseq_complete_check',
                       tags$h3('DESeq Results Exploration'),
                       uiOutput('waldContrastFactor'),
                       uiOutput('waldControlLevel'),
                       conditionalPanel('output.waldControlLevelComplete',
                                        uiOutput('waldTreatedLevel'))
                       )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel('Count and Column Data', verticalLayout(
          tags$br(),
          dataTableOutput("countDataTable"),
          tags$hr(),
          dataTableOutput("colDataTable")
        )),
        tabPanel('DESeq Analysis Results', verticalLayout(
          # verbatimTextOutput('ddsWaldPrint'),
          dataTableOutput('waldResultTable')
        )),
        tabPanel('MA plot', verticalLayout(
          tags$br(),
          # verbatimTextOutput('dfMAplot'),
          # plotlyOutput("MA_plot"),
          # tags$br(),
          # plotlyOutput("click"),
          # tags$br(),
          fluidRow(column(width=6, plotlyOutput("MA_plot")),
                   column(width=6, plotlyOutput("click"))),
          dataTableOutput("brush")
        )),
        tabPanel('Comparison'),
        tabPanel('Heatmap'),
        tabPanel('Principal Component Plot')
      )
    )
  )
)