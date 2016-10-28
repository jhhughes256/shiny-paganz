#Define server for Preliminary Patient Education Tool Application (Example 2)
#Load package libraries
	library(shiny)
	library(deSolve)
	library(ggplot2)

#Code for functions and variables which are not reactive (not dependent on "input$X")
#ggplot2 theme
	theme_custom <- theme_set(theme_grey(18))

#Function containing differential equations for amount in each compartment	
	DES <- function(T, A, THETA) {
	
		KA  <- THETA[1]
		K10 <- THETA[2]                                   			
			
        dA <- vector(length = 2)
        dA[1] = -KA*A[1]			#Absorption compartment
        dA[2] =  KA*A[1] -K10*A[2]	#Central compartment  
			
		list(dA)			
	}
		
#Make a TIME range (0 to 240 hours [10 days] at increments of 0.1 hours)
#To change the duration evaluated time-period, only the value for the "to" argument
#needs to be altered in the whole script
	TIME <- seq(from = 0, to = 240, by = 0.1)
	
#TIMElast is used in later functions for assigning dose events	
	TIMElast <- max(TIME)

#-----------------------------------------------------------------------------------
#Define user-input dependent functions for output	
	shinyServer(function(input, output) {

#Reactive expression to generate the plot
#This is called whenever the input changes
	sim.data <- reactive({
	
	#Input patient data
		AGE <- input$AGE
		WT <- input$WT
		
	#Make a parameter vector for input into DES function
	#Exponent values are large to have greater impact on the plot when input changes
		KA <- 0.5
		CL <- 10*(AGE/60)^-1.4
		V <- 175*(WT/70)^2
	
	#Calculate rate constants for the differential equations
		KA <- KA
		K10 <- CL/V
	
	#Parameter vector
		THETAlist <- c(KA,K10)
		
#----------------------------------------------------------------------------------
#Input doses specific to dosing frequency
	#Input dosing regimen
	FREQ <- input$FREQ
		
	#Input prescribed dose
	DOSE <- input$DOSE
		
#Once daily dosing		
	if (FREQ == 1) {
	
	ndoses <- TIMElast/24 + 1
		
	#Dose event data for once daily dosing (for deSolve)
	DOSEdata <- data.frame(var    = rep(1, times = ndoses),
				           time   = seq(0,TIMElast,24),
			               value  = rep(DOSE, times = ndoses),
			               method = rep("add", times = ndoses))

	#Conditions for Dose on Day 3
		#DOSE13 = Once daily dosing (1), 3rd dose (3) in 10 day period
		if (input$DOSE13 == TRUE) DOSEdata$value[DOSEdata$time == 72] <- 0
	
	#Conditions for Dose on Day 4
		#DOSE14 = Once daily dosing (1), 4th dose (4) in 10 day period
		if (input$DOSE14 == TRUE) DOSEdata$value[DOSEdata$time == 96] <- 2*DOSE
	
	}

#Twice daily dosing					
	if (FREQ == 2) {
	
	ndoses <- 2*TIMElast/24 + 1
	
	#Dose event data for twice daily dosing	
	DOSEdata <- data.frame(var    = rep(1, times = ndoses),
				           time   = seq(0,TIMElast,12),
			               value  = rep(DOSE, times = ndoses),
			               method = rep("add", times = ndoses))

	#Conditions for Second Dose on Day 2
		#DOSE26 = Twice daily dosing (2), 6th dose (6) in 10 day period
		if (input$DOSE26 == TRUE) DOSEdata$value[DOSEdata$time == 60] <- 0
		
	#Conditions for First Dose on Day 3
		#DOSE27 = Twice daily dosing (2), 7th dose (7) in 10 day period
		if (input$DOSE27 == TRUE) DOSEdata$value[DOSEdata$time == 72] <- 2*DOSE

	}
					
#Three times daily dosing					
	if (FREQ == 3) {
	
	ndoses <- 3*TIMElast/24 + 1
	
	#Dose event data for twice daily dosing	
	DOSEdata <- data.frame(var    = rep(1, times = ndoses),
				           time   = seq(0,TIMElast,8),
			               value  = rep(DOSE, times = ndoses),
			               method = rep("add", times = ndoses))

	#Conditions for Third Dose on Day 3
		#DOSE312 = Three times a day dosing (3), 12th dose (12) in 10 day period
		if (input$DOSE312 == TRUE) DOSEdata$value[DOSEdata$time == 88] <- 0
		
	#Conditions for First Dose on Day 4
		#DOSE313 = Three times a day dosing (3), 13th dose (13) in 10 day period
		if (input$DOSE313 == TRUE) DOSEdata$value[DOSEdata$time == 96] <- 2*DOSE
	
	}		
								
#Set initial conditions in each compartment
	A_0 <- c(A1 = 0, A2 = 0)
		
#Run differential equation solver (deSolve package)	
	sim.data.df <- lsoda(A_0, TIME, DES, THETAlist, events = list(data=DOSEdata))
		
#Process the simulated output	
	sim.data.df <- as.data.frame(sim.data.df)
	sim.data.df$CONC <- sim.data.df$A2/V
	sim.data.df$DAYS <- sim.data.df$time/24
		
	sim.data.df <- as.data.frame(sim.data.df)
		
	})	#Brackets closing "reactive" expression
	
#----------------------------------------------------------------------------------
#Generate a plot of the data
#Also uses the inputs to build the plot (ggplot2 package)	
	output$plotCONC <- renderPlot({	
		
		plotobj <- ggplot(sim.data())
		plotobj <- plotobj + geom_abline(aes(slope = 0, intercept = 2.5), linetype = "dashed", size = 1)
		plotobj <- plotobj + geom_abline(aes(slope = 0, intercept = 5), linetype = "dashed", size = 1)
		plotobj <- plotobj + geom_line(aes(x = DAYS, y = CONC), colour = "red", size = 1)
		plotobj <- plotobj + annotate("text", x = 9.9, y = 2.1, label = "Ineffective", colour = "black", size = 6)
		plotobj <- plotobj + annotate("text", x = 10.1, y = 5.7, label = "Toxic", colour = "black", size = 6)
		plotobj <- plotobj + scale_y_continuous("Concentration (mg/L) \n", lim = c(0,15))
		plotobj <- plotobj + scale_x_continuous("\nTime (days)", breaks = c(0,1,2,3,4,5,6,7,8,9,10))
		print(plotobj)
		
	})	#Brackets closing "renderPlot" expression
	
})	#Brackets closing "shinyServer" function
