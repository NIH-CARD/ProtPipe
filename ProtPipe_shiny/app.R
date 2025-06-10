library(shiny)

source("global.R")
source("ui.R")
source("server.R")
source("helpers.R")

shinyApp(ui = ui, server = server)
