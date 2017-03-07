library(shiny)
library(plotly)
library(DESeq2)
library(DT)



fluidPage(
  sidebarLayout(
    sidebarPanel(
      tags$h3('Data Input'),
      tags$hr(),
      fileInput('countDataInput', label='Upload count data'),
      fileInput('colDataInput', label='Upload column data')
    ),
    mainPanel(
      tabsetPanel(
        tabPanel('Count and Column Data', verticalLayout(
          tags$br(),
          dataTableOutput("countDataTable"),
          tags$hr(),
          dataTableOutput("colDataTable")
        )),
        tabPanel('MA plot', verticalLayout(
          tags$br(),
          plotlyOutput("MA_plot"),
          tags$br(),
          verbatimTextOutput("brush")
          # dataTableOutput("brush")
        ))
      )
    )
  )
)