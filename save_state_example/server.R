# server.R script for state_saving function example
# Allows user to save state of an object in the environment
# Can then be recalled or reused

# The use of "save states" or "state saving" is very useful when you want to
# reuse a previous state of your app, with the current state.
# ------------------------------------------------------------------------------
# Non-reactive objects/expressions
# Load package libraries
  library(shiny)
  library(ggplot2)
# Define custom ggplot theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
# Define sequence for time at which concentrations will be calculated
  time <- sort(unique(c(
    seq(from = 0, to = 5, by = 0.25),
    seq(from = 5, to = 24, by = 1)
  )))

# ------------------------------------------------------------------------------
# Reactive objects/expressions
  shinyServer(function(input, output) {
    # First set up data you plan to plot
    # Save data to be plotted as a reactive object due to multiple reactive
    # functions depending on it
    Rconc <- reactive({
    # Concentrations simulated from a basic one compartment model
      dose <- input$DOSE
      wt <- 70  # Weight (kg)
      crcl <- 90  # Creatinine clearance (mL/min)
      CL <- 15*((wt/70)^0.75)*((crcl/90)^1.5)  # Clearance, L/h
      V <- 20*(wt/70)  # Volume of distribution, L
      KA <- 0.2  # Absorption rate constant, h^-1
      conc <- dose*KA/V*(KA - CL/V)*(exp(-CL/V*time) - exp(-KA*time))
    })  # Rconc

    # Now set up your reactiveValue object where you will create "save states"
    # This is where the function governing the second line will be saved
    # While reactive() is used to create a reactive object, reactiveValues is
    # used to create a list of reactive values inside of a non-reactive object
    # It shares many similarities with list()
    rv <- reactiveValues(
      Sconc = NA
    )  # rv

    # The reactive function observeEvent is used to observe save button
    # On pressing the save button, save current state of plot to rv$Sconc
    observeEvent(input$save, {
      rv$Sconc <- Rconc()
    })  #observeEvent

    # Observe the clear button for when the user wishes to delete the saved plot
    observeEvent(input$clear, {
      rv$Sconc <- NA
    })  #observeEvent

    # There are now two concentration vectors
    #   Rconc - reactive concentrations - called with Rconc()
    #   Sconc - saved concentrations - called with rv$Sconc

    # The server is now ready to use "save states"
    # Render your plot (table/text/etc.)
    output$concPlot <- renderPlot({
    # Create the basic ggplot2 object for output
      plotobj <- ggplot()
      plotobj <- plotobj + geom_line(aes(x = time, y = Rconc()), colour = "red")
      plotobj <- plotobj + scale_x_continuous("\nTime (hours)")
      plotobj <- plotobj + scale_y_continuous("Concentration (mg/L)\n",
        lim = c(0,0.3))
    # If statement to also plot the "save state" if it exists
    # Index for first number used to silence errors inherent with using is.na
    # with vector length more than 1
      if (!is.na(rv$Sconc[1])) {
        plotobj <- plotobj + geom_line(aes(x = time, y = rv$Sconc),
          colour = "blue")
      }
    # Print final object
      plotobj
    })  # renderPlot
  })  # shinyServer
