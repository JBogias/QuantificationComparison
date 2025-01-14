library(shiny)
library(here)

ui <- shinyUI(
  fluidPage(
    includeHTML(here("www/filtering.html"))
  )
)

server <- function(input, output) {}

shinyApp(ui, server)

