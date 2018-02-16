#====== nested R functions that assemble HTML user interface ========

fluidPage(
  sidebarLayout(
    sidebarPanel(
      titlePanel('Input Types'),
      numericInput(inputId = 'numericInput', label = "1. numeric input", value = 25),
      actionButton(inputId = 'actionButton',label = '2. action button'),
      tags$br(),
      tags$br(),
      actionLink(inputId = 'actionLink', label = '3. action Link'),
      checkboxGroupInput(inputId = 'checkboxGroupInput', label = '4. checkbox group input', choices = c('a for apple', 'b for ball'), selected = ('b for ball')),
      dateInput(inputId = 'dateInput', label = '5. date input'),
      dateRangeInput(inputId = 'dateRangeInput', label = '6. date range input'),
      fileInput(inputId = 'fileInput', label = '7. file input'),
      passwordInput(inputId = 'passwordInput', label = '8. password input')
    ),
    mainPanel(
      titlePanel('Input Values'),
      verbatimTextOutput(outputId = "o_numericInput"),
      verbatimTextOutput(outputId = "o_actionButton"),
      verbatimTextOutput(outputId = "o_actionLink"),
      verbatimTextOutput(outputId = "o_checkoutGroupInpu"),
      verbatimTextOutput(outputId = "o_dateInput"),
      verbatimTextOutput(outputId = "o_dateRangeInput"),
      verbatimTextOutput(outputId = 'o_fileInput'),
      verbatimTextOutput(outputId = 'o_passwordInput')
    )
  )
)