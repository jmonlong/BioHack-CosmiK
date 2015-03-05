library(shiny)
library(animint)

shinyUI(fluidPage(
  headerPanel("Interactive clusters"),
  sidebarPanel(
    sliderInput("nb.units", "Number of units", 10, 30, 10, step = 5),
    radioButtons("col", "Color", c("clust","dist.neighbours","nb.genes", "med.exp", "sd.exp"), selected = "nb.genes")),
 ## Main plot
  mainPanel(
    tabsetPanel(
      tabPanel("SOM", plotOutput("SOM")),
      tabPanel("Animation", animintOutput("anim"))
    )
  )))
