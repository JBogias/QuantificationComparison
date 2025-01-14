library(shiny)
library(here)

ui <- shinyUI(
  fluidPage(
    includeHTML(here("www/kmer_analysis.html"))
  )
)

server <- function(input, output) {}

shinyApp(ui, server)