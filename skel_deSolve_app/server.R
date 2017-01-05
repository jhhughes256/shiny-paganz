# R script for simulating a population and concentrations as described by:
# 2-compartment model
# Dosing regimen:
#   IV bolus    (at 0 hours)
#   Oral Dosing (at 1 hours, then every 12 hours after IV bolus)
#   IV Infusion (at 72 hours, duration 10 hours)
# Population parameter variability on CL, V1, KA
# Variance-covariance matrix for CL, V1, KA
# Covariate effects:
#   Gender & Creatinine Clearance on CL

# Proportional error model (with optional additive residual)
# ------------------------------------------------------------------------------
# Define user-input dependent functions for output
shinyServer(function(input, output) {
# Reactive expression to generate a reactive data frame
# Whenever an input changes this function re-evaluates
	Rpar.data <- reactive({
	# Create a parameter dataframe with ID and parameter values for each individual
	# Define individual
	  n <- 10  # Number of "individuals"
	  ID <- seq(from = 1, to = n, by = 1)  # Simulation ID
	  WT <- input$wt  # Total body weight, kg
	  AGE <- input$age  # Age, years
	  SECR <- 60  # Serum Creatinine, umol/L
	  SEX <- 1  # Gender, Male (0) Female (1)
	  SMOK <- 0  # Smoking Status, Not Current (0) Current (1)

	# Now use multivariate rnorm to turn the covariance matrix into ETA values
	# Generate random samples from multivariate distribution
	  ETAmat <- mvrnorm(n = n, mu = c(0, 0, 0), OMEGA)
	  ETA1 <- ETAmat[, 1]  # Take first column of matrix
	  ETA2 <- ETAmat[, 2]
	  ETA3 <- ETAmat[, 3]

	# Define covariate effects
	  SMOKCOV <- 1
	  if (SMOK == 1) SMOKCOV <- SMOKCOV + COV1
	  CRCL <- ((140 - AGE)*WT)/(SECR*0.815)  # Male creatinine clearance
	  if (SEX == 0) CRCL <- CRCL*0.85  # Female creatinine clearance

	# Define individual parameter values
	  CL <- CLPOP*exp(ETA1)*((WT/70)^0.75)*SMOKCOV*((CRCL/90)^COV2)
	  V1 <- V1POP*exp(ETA2)*(WT/70)
	  Q  <- QPOP*(WT/70)^0.75
	  V2 <- V2POP*(WT/70)
	  KA <- KAPOP*exp(ETA3)

	# Calculate rate-constants for differential equation solver
	  K12 <- Q/V1
	  K21 <- Q/V2
	  K10 <- CL/V1

	# Collect the individual parameter values in a parameter dataframe
	  par.data <- data.frame(
	    ID, CL, V1, Q, V2,  # Patient parameters
	    KA, K12, K21, K10,  # Rate constants
	    WT, AGE, SECR, SEX, SMOK # Covariates
		)
	})  # Rpar.data

#------------------------------------------------------------------------------
	Rsim.data <- reactive({

		##########
		##_ORAL_##
		##########

		########
		##_IV_##
		########

	# Set blank event data for deSolve
	# Currently no IV bolus or oral dosing events
		all.dose.data <- data.frame(var = 1, time = 0, value = 0, method = "add")

	# Define infusion using input parameters
	# If statement for checkbox ui
		if (input$infcheck == FALSE) {
			inf.dose <- 0
			inf.dur <- 0
			inf.start <- 0
		} else {
			inf.dose <- input$infdose  #mg
		  inf.dur <- as.numeric(input$infdur)
			inf.start <- input$inftimes
		}
	  inf.rate <- inf.dose/inf.dur
		inf.times <- c(inf.start, inf.start + inf.dur)

	# Make a time sequence (hours)
		all.times <- sort(unique(c(set.time, inf.times))) #, oral.dose.times, iv.dose.times)))
		final.time <- max(all.times)
		# The time sequence must include all "event" times for deSolve, so added here
		# Do not repeat a time so use "unique" as well

	# Specify infusion event times and corresponding rates
	# Note - inf.times vector of 2 values, inf.rate vector of 1 value
    inf.time.data <- c(0, inf.times, final.time)
		inf.rate.data <- c(0, inf.rate, 0, 0)

#------------------------------------------------------------------------------
	# Apply simulate.conc.cmpf to each individual in par.data
	# Maintain their individual values for V1 for later calculations
	  sim.data <- ddply(Rpar.data(), .(ID, V1), simulate.conc.cmpf,
			times = all.times, event.data = all.dose.data,
			inf.rate.fun = inf.rate.fun(inf.time.data, inf.rate.data)
		)

	# Calculate individual concentrations in the central compartment
	  sim.data$IPRED <- sim.data$A2/sim.data$V1
		return(sim.data)
	})	#Rsim.data

#-------------------------------------------------------------------------------
# Generate a plot of the data
# Also uses the inputs to build the plot
	output$plotconc <- renderPlot({
	# Generate a plot of the sim.data
	  plotobj <- NULL
	  plotobj <- ggplot(data = Rsim.data())

		########
		##_CI_##
		########
	# Shaded ribbon for prediction intervals
		plotobj <- plotobj + stat_summary(aes(x = time, y = IPRED),
			fun.ymin = CI80lo, fun.ymax = CI80hi, geom = "ribbon",
			fill = "red", alpha = 0.3)

	# Solid line for median
	  plotobj <- plotobj + stat_summary(aes(x = time, y = IPRED),
	    fun.y = median, geom = "line", size = 1, colour = "red")

		#########
		##_LOG_##
		#########

	# Continuous y axis
		plotobj <- plotobj + scale_y_continuous("Concentration (mg/L) \n",
    breaks = seq(from = 0, to = max(Rsim.data()$IPRED),
			by = ceiling(max(Rsim.data()$IPRED)/10)), lim = c(0, max(Rsim.data()$IPRED)))

	# Continuous x axis
	  plotobj <- plotobj + scale_x_continuous("\nTime (hours)", lim = c(0, 120),
			breaks = seq(from = 0,to = 120,by = 24))

	  print(plotobj)
	})	#renderPlot

# Render text for infusion rate calculation
	output$infrate <- renderText({
		inf.dur <- as.numeric(input$infdur)
		paste0("Infusion rate = ", signif(input$infdose/inf.dur, digits = 3) ," mg/hr")
	})	#renderText

	##########
  ##_CRCL_##
  ##########

})	#shinyServer
