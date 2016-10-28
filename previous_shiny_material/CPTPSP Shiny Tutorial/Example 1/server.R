#Define server for Flip-Flop Kinetics Application
#Load shiny and ggplot2 packages
	library(shiny)
	library(ggplot2)

#Code for functions and variables which are not reactive (not dependent on "input$X")	
#ggplot2 theme	
	theme_custom <- theme_set(theme_grey(18))

#TIME sequence for concentrations to be calculated
	TIME <- seq(from = 0, to = 24, by = 0.1)

#Define dose (D)
	D <- 500
		
#----------------------------------------------------------------------------------------
#Define user-input dependent functions for output
	shinyServer(function(input, output) {

#Generate a plot of the data
#Also uses the inputs to build the plot (ggplot2 package)
	output$plotCONC <- renderPlot({
	
		#Collect input from user-widgets
		KA <- input$KA
		KE <- input$KE
		V  <- input$V		
		
		#Calculate function of Dose, V, KA and KE (where F = 1)
		if (KA==KE) A <- 0 else A <- (D*KA)/(V*(KA-KE))

		#Calculate concentration
		CONC <- A*exp(-KE*TIME) - A*exp(-KA*TIME)

		#Create a dataframe of time and concentration data
		conc.data <- data.frame(TIME = TIME, CONC = CONC)
		
		#ggplot2 object
		plotobj <- ggplot(data = conc.data)
		plotobj <- plotobj + geom_line(aes(x = TIME, y = CONC), colour = "red", size = 1)
		plotobj <- plotobj + scale_y_continuous("Concentration (mg/L) \n", lim = c(0,50))
		plotobj <- plotobj + scale_x_continuous("\nTime (hours)", lim = c(0,24))
		print(plotobj)
		
	})	#Brackets closing "renderPlot"
	
})	#Brackets closing "shinyServer"
