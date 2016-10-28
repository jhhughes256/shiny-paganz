#ui.R script for shinyhangout_example appplication
#Defines widgets for reactive input for server.R and calls output objects from server.R and presents them in the user-interface
#----------------------------------------------------------------------------------
#All user-interface elements need to be bound within a page formatting function, i.e., "fixedPage"
#This application places each element below the previous one in order from top to bottom
fixedPage(
  h3("1-compartment, first-order oral absorption kinetics"),  #Application title

  hr(),  #Horizontal separating line

  fixedRow(column(6,
                  plotOutput("concPlot",width = 600)  #Concentration-time profile output
  ),
  column(6,  selectInput("DOSE",
                         "Dose (mg):",
                         choices = list("50 mg" = 50,"100 mg" = 100),
                         selected = 1),  #Select input for dose

         sliderInput("WT",
                     "Weight (kg):",
                     min = 40,
                     max = 100,
                     step = 1,
                     value = 70))
                  ),

  align = "center"  #Align all elements in the centre of the application

)  #Brackets closing "fixedPage"
