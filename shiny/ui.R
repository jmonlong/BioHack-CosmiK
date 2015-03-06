library(shiny)
library(animint)

shinyUI(fluidPage(
  headerPanel("Interactive clusters"),
  sidebarPanel(
    numericInput("nb.samp", "Number of samples", 10, 10, 30, step = 5),
    numericInput("nb.iter", "Number of iterations", 10, 10, 100, step = 10),
    numericInput("nb.units", "Number of units", 10, 10, 40, step = 5),
    numericInput("nb.clust", "Number of clusters", 2, 2,20, step = 1),
    radioButtons("med.norm", "Median normalization per gene", c(FALSE,TRUE), selected = FALSE),
    selectInput("col","Color by", c(cluster="clust",distance="dist.neighbours",genes="nb.genes", expression="med.exp", variation="sd.exp"), selected = "nb.genes")),
 ## Main plot
  mainPanel(
    tabsetPanel(
      tabPanel("SOM", column(5,plotOutput("som")), column(5, plotOutput("som.conv"))),
      tabPanel("SOM Weights", plotOutput("som.codes", height="800px")),
      tabPanel("Gene exploration", animintOutput("anim"))
    )
  )))
