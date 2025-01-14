library(shiny)
library(here)

ui <- shinyUI(
  fluidPage(
    includeHTML(here("www/annotations.html"))
  )
)

server <- function(input, output) {}

shinyApp(ui, server)