library(shiny)
fileUploadUI <- function(id, label = "Upload File") {
  ns <- NS(id)
  tagList(
    fileInput(ns("file"), label),
    actionButton(ns("clear"), "Remove file")
  )
}

fileUploadServer <- function(id, label = "Upload File") {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    file <- reactiveVal(NULL)

    output$file_ui <- renderUI({
      fileInput(ns("file"), label)
    })

    observeEvent(input$file, {
      file(input$file)
    })

    observeEvent(input$clear, {
      file(NULL)
      # Force UI to re-render, fully resetting the file input
      output$file_ui <- renderUI({
        fileInput(ns("file"), label)
      })
    })

    return(file)
  })
}

