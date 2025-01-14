library(shiny)
library(here)

ui <- shinyUI(
  fluidPage(
    includeHTML(here("www/transposons.html"))
  )
)

server <- function(input, output) {}

shinyApp(ui, server)