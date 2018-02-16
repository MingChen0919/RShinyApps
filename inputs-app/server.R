#===== a function with instructions on how to build and rebuild the R objects ======

function(input, output) {
  
  output$o_numericInput = renderPrint({cat('1. numericInput \n'); str(input$numericInput)})
  output$o_actionButton = renderPrint({cat('2. actionButton \n'); str(input$actionButton)})   
  output$o_actionLink = renderPrint({cat('3. actionLink \n'); str(input$actionLink)})
  output$o_checkoutGroupInpu = renderPrint({cat('4. checkoutGroupInput \n'); str(input$checkboxGroupInput)})
  output$o_dateInput = renderPrint({cat('5. dateInput \n'); str(input$dateInput)})
  output$o_dateRangeInput = renderPrint({cat('6. dateRangeInput \n'); str(input$dateRangeInput)})
  output$o_fileInput = renderPrint({cat('7. fileInput \n'); str(input$fileInput)})
  output$o_passwordInput = renderPrint({cat('8. passwordInput \n'); str(input$passwordInput)})
  output$o_radioButtons = renderPrint({cat('9. radioButtons \n'); str(input$radioButtons)})
  output$o_selectInput = renderPrint({cat('10. selectInput \n'); str(input$selectInput)})
  output$o_sliderInput = renderPrint({cat('11. sliderInput \n'); str(input$sliderInput)})
  output$o_textInput = renderPrint({cat('12. textInput \n'); str(input$textInput)})
  # output$o_submitButton = renderPrint({})
}