#ui.R script for shinyhangout_example appplication
#Defines widgets for reactive input for server.R and calls output objects from server.R and presents them in the user-interface
#----------------------------------------------------------------------------------
#All user-interface elements need to be bound within a page formatting function, i.e., "fixedPage"
#This application places each element below the previous one in order from top to bottom
fixedPage(
  h3("1-compartment, first-order oral absorption kinetics"),  #Application title

  hr(),  #Horizontal separating line

  plotOutput("concPlot",width = 600),  #Concentration-time profile output

  sliderInput("WT",
              "Weight (kg):",
              min = 40,
              max = 100,
              step = 5,
              value = 70),  #Slider input for weight

  sliderInput("CRCL",
              "Creatinine clearance (mL/min):",
              min = 30,
              max = 90,
              step = 10,
              value = 60),  #Slider input for creatinine clearance

  selectInput("DOSE",
              "Dose administered (mg):",
              choices = list("50 mg" = 1,"100 mg" = 2),
              selected = 1),  #Select box input for dose

  checkboxInput("LOGSCALE",
                "Plot concentrations on a log-scale:",
                value = FALSE),  #Checkbox input for y-axis scale

  align = "center"  #Align all elements in the centre of the application

)  #Brackets closing "fixedPage"
