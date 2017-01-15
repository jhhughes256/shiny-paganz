# PAGANZ 2017 - Introduction to the Shiny Framework for Pharmacometricians - Theory
# ui.R script for theory_basic_example application
# Defines widgets for reactive input for server.R and calls output objects from server.R and presents them in the user-interface
# ----------------------------------------------------------------------------------
# All user-interface elements need to be bound within a page formatting function, i.e., "fixedPage"
# This application places each element below the previous one in order from top to bottom
  fixedPage(
    h3("1-compartment, first-order oral absorption kinetics"),  # Application title

    hr(),  # Horizontal separating line

    plotOutput("concPlot",width = 600),  # Concentration-time profile output

    sliderInput("DOSE",
      "Dose (mg):",
      min = 10,
      max = 100,
      step = 5,
      value = 50
    ), # Slider input for dose

    align = "center"  # Align all elements in the centre of the application

  )  # fixedPage
