library(shiny)

shinyUI(fluidPage(
  headerPanel("Interactive clusters"),
  sidebarPanel(
    shiny::tabPanel("compile", shiny::actionButton("compile","Compile animation"), textOutput("compile")),
    sliderInput("nb.units", "Number of units", 10, 30, 10, step = 5),
    radioButtons("col", "Color", c("clust","dist.neighbours","nb.genes", "med.exp", "sd.exp"), selected = "nb.genes")),
 ## Main plot
  mainPanel(
    tags$head(includeScript('www/js/animint.js')),
    tags$head(includeScript('www/js/d3.v3.js')),
    shiny::tabsetPanel(
      shiny::tabPanel("SOM", shiny::plotOutput("SOM")),
      shiny::tabPanel("Animation", includeHTML("index.html"))
    )
  )))
