library(shiny)
library(here)

ui <- shinyUI(
  fluidPage(
    includeHTML(here("www/create_dtelists.html"))
  )
)

server <- function(input, output) {}

shinyApp(ui, server)
