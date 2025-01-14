# App for ACR and App Service deployment
# I would just use the docker-compose for ACR but then I'd need to store 
# all seven of the images in the ACR, which would rake up costs
# App service is free though
library(shiny)
library(here)

addResourcePath(prefix = "html", directoryPath = "www/")

ui <- fluidPage(
  class = "fluid-page",
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "html/styling.css")
  ),
  titlePanel(paste0("Code Companion for my Multi-Mapping Analysis workflow",
                    " by Justin Bogias")),
  tabsetPanel(
    tabPanel("Prepare Annotations",
             class = "tab-panel",
             tags$iframe(src = "html/annotations.html",
                         style = "width:70vw;height:80vh;",
                         class = "markdown_frame")),
    tabPanel("Prepare DTELists",
             class = "tab-panel",
             tags$iframe(src = "html/create_dtelists.html",
                         style = "width:70vw;height:80vh;",
                         class = "markdown_frame")),
    tabPanel("Filtering",
             class = "tab-panel",
             tags$iframe(
               src = "html/filtering.html",
               style = "width:70vw;height:80vh;",
               class = "markdown_frame")),
    tabPanel("Principal Component Analysis",
             class = "tab-panel",
             tags$iframe(src = "html/principal_component_analysis.html",
                         style = "width:70vw;height:80vh;",
                         class = "markdown_frame")),
    tabPanel("DTE and DTU",
             class = "tab-panel",
             tags$iframe(src = "html/DTE_DTU.html",
                         style = "width:70vw;height:80vh;",
                         class = "markdown_frame")),
    tabPanel("Correlation",
             class = "tab-panel",
             tags$iframe(
               src = "html/06_correlations.html",
               style = "width:70vw;height:80vh;",
               class = "markdown_frame")),
    tabPanel("Multi-Mapping Analysis",
             class = "tab-panel",
             tags$iframe(
               src = "html/multimap_analysis.html",
               style = "width:70vw;height:80vh;",
               class = "markdown_frame")),
    tabPanel("K-mer Analysis",
             class = "tab-panel",
             tags$iframe(
               src = "html/kmer_analysis.html",
               style = "width:70vw;height:80vh;",
               class = "markdown_frame")),
    tabPanel("Transposons",
             class = "tab-panel",
             tags$iframe(
               src = "html/transposons.html",
               style = "width:70vw;height:80vh;",
               class = "markdown_frame")),
    tabPanel("Plot Annotations",
             class = "tab-panel",
             tags$iframe(
               src = "html/plot_annotations.html",
               style = "width:70vw;height:80vh;",
               class = "markdown_frame")),
    tabPanel("Plot PCBP2",
             class = "tab-panel",
             tags$iframe(
               src = "html/pcbp2_plotting.html",
               style = "width:70vw;height:80vh;",
               class = "markdown_frame"))
  )
)

server <- function(input, output, session) {
}

shinyApp(ui, server)