library(shiny)
library(expm)
library(plotly)
library(ggplot2)
library(ggplotify)
library(graphics)


# logifySlider javascript function
JS.logify <-
"
// function to logify a sliderInput
function logifySlider (sliderId, sci = false) {
  if (sci) {
    // scientific style
    $('#'+sliderId).data('ionRangeSlider').update({
      'prettify': function (num) { return ('10<sup>'+num+'</sup>'); }
    })
  } else {
    // regular number style
    $('#'+sliderId).data('ionRangeSlider').update({
      'prettify': function (num) { return (Math.pow(10, num)); }
    })
  }
}"

# call logifySlider for each relevant sliderInput
JS.onload <-
"
// execute upon document loading
$(document).ready(function() {
  // wait a few ms to allow other scripts to execute
  setTimeout(function() {
    // include call for each slider
    logifySlider('nsites', sci = false)
    logifySlider('popsize', sci = false)
    logifySlider('log_slider2', sci = true)
  }, 5)})
"


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  tags$head(tags$script(HTML(JS.logify))),
  tags$head(tags$script(HTML(JS.onload))),

  # App title ----
  titlePanel("Probability of having a recombination event within a locus"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Slider for the number of sites ----
      sliderInput(inputId = "nsites",
                  label = "Number of sites:",
                  min = 0,
                  max = 6,
                  value = 3,
                  step = 1),

      # Input: Slider for the number of sites ----
      sliderInput(inputId = "popsize",
                  label = "Population size:",
                  min = 0,
                  max = 6,
                  value = 3,
                  step = 1)
    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Histogram ----
      plotlyOutput(outputId = "distPlot")

    )
  )
)

prob <- function(t,rho,theta) {
    theta / (theta+rho) + rho / (theta+rho)*exp(-(theta+rho)*t)
}

# Define server logic required to draw a histogram ----
server <- function(input, output) {

  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$distPlot <- renderPlotly({


  # theta    = 1.0 / pop_size / 2.0;
  # rho      = RECOMB_RATE * genome_size;
  # analytical_prob = theta / (rho+theta) + rho / (rho+theta) * std::exp(-(rho+theta)*time);

  generations = 10^(0:6)
  pop_size    = 10^input$popsize
  num_sites   = 10^input$nsites
  RECOMB_RATE = 2.52 * 1E-8

  rho      = RECOMB_RATE * num_sites
  theta    = 1.0 / pop_size / 2.0
  time    = generations
  # analytical solution in the 2 sites scenario
  #  p <- 1 - (theta / (rho+theta) + rho / (rho+theta) * exp(-(rho+theta)*time))
  max_dim <- 10
  x <- matrix(rep(0,max_dim*max_dim), max_dim, max_dim)
  for (k in 1:max_dim) {
    if ( k < max_dim ) {
      j = k+1
      x[k,j] <- (num_sites-k)*RECOMB_RATE
      x[k,k] <- x[k,k] - (num_sites-k)*RECOMB_RATE
    }

    if ( k > 1 ) {
      j = k-1
      x[k,j] <- k*(k-1.0)/2.0/pop_size/2.0
      x[k,k] <- x[k,k] - k*(k-1.0)/2.0/pop_size/2.0
    }
  }

  yyy <- c()
  for ( k in 1:length(generations) ) {
    tp     <- expm(x*generations[k])
    yyy[k] <- 1 - tp[1,1]
  }
  
  df = data.frame(x=generations, y=yyy)
  p <- ggplot(df,aes(x=x,y=y))+geom_point()+
    scale_x_log10()+ylim(c(0,1))
  
  ggplotly(p)

  })

}



shinyApp(ui = ui, server = server)
