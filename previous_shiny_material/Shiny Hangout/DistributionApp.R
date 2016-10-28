#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)

theme <- theme_set(theme_bw(base_size = 14))

# Define UI for application that draws a inverse gamma distribution
ui <- shinyUI(fluidPage(

   # Application title
   titlePanel("Distributions"),

   # Sidebar with a slider input for gamma parameters and length of x
   sidebarLayout(
      sidebarPanel(
        selectInput("disttype","Distribution type",choices = list("Normal" = 1,"Log-Normal" = 2,"Inverse Gamma" = 3,"Cauchy" = 4,"Uniform" = 5)),
        sliderInput("lengthx","Length of x:",min = 0,max = 10,step = 1,value = 10),
        conditionalPanel(condition = "input.disttype == 1",
          sliderInput("nmean","Mean:",min = -10,max = 10,step = 0.01,value = 0)
        ),
        conditionalPanel(condition = "input.disttype == 2",
          sliderInput("lnmean","Mean:",min = 0.01,max = 20,value = 1)
        ),
        conditionalPanel(condition = "input.disttype <= 2",
          sliderInput("sd","Standard Deviation:",min = 0,max = 10,step = 0.01,value = 1)
        ),
        conditionalPanel(condition = "input.disttype == 3",
          sliderInput("alpha","alpha:",min = 0.01,max = 2,step = 0.01,value = 0.01),
          sliderInput("beta","beta:",min = 0.01,max = 2,step = 0.01,value = 0.01)
        ),
        conditionalPanel(condition = "input.disttype == 4",
          sliderInput("l","Location:",min = 0.01,max = 2,step = 0.01,value = 0.01),
          sliderInput("s","Scale:",min = 0.01,max = 2,step = 0.01,value = 0.01)
        ),
        conditionalPanel(condition = "input.disttype == 5",
          sliderInput("min","Minimum:",min = 0,max = 10,step = 0.5,value = 1),
          sliderInput("max","Maximum:",min = 0,max = 10,step = 0.5,value = 7)
        )
      ),

      # Show a plot of the generated distribution
      mainPanel(
        plotOutput("distPlot")
      )
   )
))

# Define server logic required to the inverse gamma distribution
server <- shinyServer(function(input, output) {

   output$distPlot <- renderPlot({
      #Generate a sequence of x based on "lengthx" from the ui.R
      if (input$disttype == 1) {
        x <- seq(from = -10,to = input$lengthx,by = 0.001)
      }
      else if (input$disttype == 2) {
        x <- seq(from = 0.001,to = input$lengthx,by = 0.001)
      } else {
        x <- seq(from = 0,to = input$lengthx,by = 0.001)
      }
      #Generate the normal distribution for each value of x given the input for mean and standard deviation
      if (input$disttype == 1) {
        dist <- 1/(sqrt(2*pi)*input$sd)*exp(-((x-input$nmean)^2)/(2*input$sd^2))
      }
      #Generate teh log-normal distribution for each value of x given the input for mean and standard deviation
      if (input$disttype == 2) {
        dist <- 1/(sqrt(2*pi)*input$sd)*exp(-((log(x)-log(input$lnmean))^2)/(2*input$sd^2))
      }
      #Generate the inverse gamma distribution for each value of x given the input for alpha and beta
      if (input$disttype == 3) {
        dist <- (input$beta^input$alpha)/gamma(input$alpha)*x^(-input$alpha-1)*exp(-input$beta/x)
      }
      #Generate the cauchy distribution for each value of x given the input for l (location) and s (scale)
      if (input$disttype == 4) {
        dist <- 1/(pi*input$s*(1+((x-input$l)/input$s)^2))
      }
      #Generate the uniform distribution for each value of x given the input for min and max
      if (input$disttype == 5) {
        dist <- x/x*1/(input$max-input$min)
      }
      #Plot the distribution
      plotobj1 <- NULL
      plotobj1 <- ggplot()
      plotobj1 <- plotobj1 + geom_line(aes(x = x,y = dist))
      plotobj1 <- plotobj1 + scale_y_continuous("Distribution Density\n")
      plotobj1 <- plotobj1 + scale_x_continuous("\nx")
      print(plotobj1)
   })
})

# Run the application
shinyApp(ui = ui, server = server)

