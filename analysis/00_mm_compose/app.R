library(shiny)

ui <- fluidPage(
  titlePanel("Gestational Diabetes Mellitus"),
  tabsetPanel(
    tabPanel("Load Annotations",
             tags$iframe(src = "http://localhost:4001",
                         style = "width:70vw;height:100vh;")),
    tabPanel("Prepare DTEList", 
             tags$iframe(src = "http://localhost:4002",
                         style = "width:70vw;height:100vh;")),
    tabPanel("Filtering",
             tags$iframe(src = "http://localhost:4003",
                         style = "width:70vw;height:100vh;")),
    tabPanel("Principal Component Analysis",
             tags$iframe(src = "http://localhost:4004",
                         style = "width:70vw;height:100vh;")),
    tabPanel("DTE and DTU Analysis",
             tags$iframe(src = "http://localhost:4005",
                         style = "width:70vw;height:100vh;")),
    tabPanel("Correlations",
             tags$iframe(src = "http://localhost:4006",
                         style = "width:70vw;height:100vh;")),
    tabPanel("Multi-Mapping Analysis",
             tags$iframe(src = "http://localhost:4007",
                         style = "width:70vw;height:100vh;")),
    tabPanel("K-mer Analysis",
             tags$iframe(src = "http://localhost:4008",
                         style = "width:70vw;height:100vh;")),
    tabPanel("Transposon Analysis", 
             tags$iframe(src = "http://localhost:4009",
                         style = "width:70vw;height:100vh;")),
    tabPanel("Plot Annotations",
             tags$iframe(src = "http://localhost:4010",
                         style = "width:70vw;height:100vh;")),
    tabPanel("PCBP2 Plotting",
             tags$iframe(src = "http://localhost:4011",
                         style = "width:70vw;height:100vh;"))
  )
)

server <- function(input, output, session) {
}

shinyApp(ui, server)