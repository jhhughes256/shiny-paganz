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
TIME <- sort( unique(c(seq(from = 0,to = 5,by = 0.25),seq(from = 5,to = 24,by = 1))))

#----------------------------------------------------------------------------------
#Reactive objects/expressions, i.e., those that depend on values from the ui.R need to be within "shinyServer"
shinyServer(function(input, output) {  #Requires "input" and "output" objects
  output$concPlot <- renderPlot({  #"renderPlot" creates a "reactive" plot output object called "concPlot"

    #Assign the stored widget values from ui.R to objects in the "renderPlot" expression
    WT <- input$WT  #Call in the widget value for "WT", weight (kg)
    CRCL <- input$CRCL  #Call in the widget value for "CRCL", creatinine clearance (mL/min)
    if (input$DOSE == 1) {  #Call in the widget value for "DOSE"
      DOSE <- 50  #Where option "1" was the 50 mg dose
    } else {
      DOSE <- 100  #Where option "2" was the 100 mg dose
    }

    #Calculate pharmacokinetic parameters based on widget input and model's structure
    CL <- 15*(WT/70)^0.75*(CRCL/90)^-1.5  #Clearance, L/h, dependent on weight and creatinine clearance input
    V <- 20*(WT/70)  #Volume of distribution, L, dependent on weight input
    KA <- 0.2  #Absorption rate constant, h^-1
    #For each value of TIME (described earlier), simulate concentrations based on the above pharmacokinetic parameters
    #1-compartment, first-order oral absorption model
    CONC <- DOSE*KA/V*(KA-CL/V)*(exp(-CL/V*TIME)-exp(-KA*TIME))

    #Create the ggplot2 object for output
    plotobj <- ggplot()
    plotobj <- plotobj + geom_line(aes(x = TIME,y = CONC),colour = "red")
    plotobj <- plotobj + scale_x_continuous("\nTime (hours)")
    if (input$LOGSCALE == FALSE) {  #If the "LOGSCALE" checkbox is unticked, then plot on a linear-scale
      plotobj <- plotobj + scale_y_continuous("Concentration (mg/L)\n",lim = c(0,5))
    }
    if (input$LOGSCALE == TRUE) {  #If the "LOGSCALE" checkbox is ticked, then plot on a log-scale
      plotobj <- plotobj + scale_y_log10("Concentration (mg/L)\n",lim = c(0.001,5))
    }
    print(plotobj)  #This is the resulting object that will be sent to ui.R

  })  #Brackets closing "renderPlot" expression
})  #Brackets closing "shinyServer" function
