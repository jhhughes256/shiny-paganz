# ui.R script for state_saving function example
# ------------------------------------------------------------------------------
  fixedPage(
    h3("1-compartment, first-order oral absorption kinetics"),
    hr(),
    plotOutput("concPlot",width = 600),  # Concentration-time profile output
    column(9,
      sliderInput("DOSE",
        "Dose (mg):",
        min = 10,
        max = 100,
        step = 5,
        value = 50
      )  # sliderInput
    ),  # column
    column(2,
      actionButton("save", "Save"),  # actionButton
      actionButton("clear", "Clear")  # actionButton
    ),  # column
    align = "center"
  )  # fixedPage
