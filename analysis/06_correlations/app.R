library(shiny)
library(here)

ui <- shinyUI(
  fluidPage(
    includeHTML(here("www/06_correlations.html"))
  )
)

server <- function(input, output) {}

shinyApp(ui, server)