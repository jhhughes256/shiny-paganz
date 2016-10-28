#server.R script for shinyhangout_example application
#Calls reactive input from ui.R to simulate a concentration-time profile and plots the resulting profile
#----------------------------------------------------------------------------------
#Non-reactive objects/expressions, i.e., do not depend on reactive widget input from ui.R
#Load package libraries
library(shiny)  #Package for creating user-interface for R program
library(ggplot2)  #Plotting

#Define custom ggplot theme
theme_bw2 <- theme_set(theme_bw(base_size = 14))

#Define sequence for time at which concentrations will be calculated
#Heavier sampling at the beginning of the sequence (every 15 minutes), then every hour after 5 hours since administration
TIME <- sort(unique(c(seq(from = 0,to = 5,by = 0.25),seq(from = 5,to = 24,by = 1))))

#----------------------------------------------------------------------------------
#Reactive objects/expressions, i.e., those that depend on values from the ui.R need to be within "shinyServer"
shinyServer(function(input, output) {  #Requires "input" and "output" objects
  output$concPlot <- renderPlot({  #"renderPlot" creates a "reactive" plot output object called "concPlot"

    #Assign the stored widget values from ui.R to objects in the "renderPlot" expression
    DOSE <- as.numeric(input$DOSE)
    #Other covariate information
    WT <- input$WT  #Weight (kg)
    CRCL <- 90  #Creatinine clearance (mL/min)

    #Calculate pharmacokinetic parameters based on widget input and model's structure
    CL <- 15*((WT/70)^0.75)*((CRCL/90)^1.5)  #Clearance, L/h
    V <- 20*(WT/70)  #Volume of distribution, L
    KA <- 0.2  #Absorption rate constant, h^-1
    #For each value of TIME (described earlier), simulate concentrations based on the above pharmacokinetic parameters
    #1-compartment, first-order oral absorption model
    CONC <- DOSE*KA/V*(KA-CL/V)*(exp(-CL/V*TIME)-exp(-KA*TIME))

    #Create the ggplot2 object for output
    plotobj <- ggplot()
    plotobj <- plotobj + geom_line(aes(x = TIME,y = CONC),colour = "red")
    plotobj <- plotobj + scale_x_continuous("\nTime (hours)")
    plotobj <- plotobj + scale_y_continuous("Concentration (mg/L)\n",lim = c(0,0.3))
    print(plotobj)  #This is the resulting object that will be sent to ui.R
    print(head(CONC))
    print(plotobj)

  })  #Brackets closing "renderPlot" expression
})  #Brackets closing "shinyServer" function
