library(shiny)
library(here)

ui <- shinyUI(
  fluidPage(
    includeHTML(here("www/principal_component_analysis.html"))
  )
)

server <- function(input, output) {}

shinyApp(ui, server)

