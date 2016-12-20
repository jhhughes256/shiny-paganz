# Globally used objects and functions
# ------------------------------------------------------------------------------
# Load package libraries
	library(shiny)	# Interactive applications
	library(ggplot2)	# Plotting
	library(plyr)  # Apply a function to a data.frame
	library(deSolve)	# Differential equation solver
	library(MASS)  # Used for mvrnorm (random samples from multivariate distribution)
	library(MBESS)  # Used for cor2cov (convert from correlation to covariance matrix)
	library(compiler)	# Compile repeatedly used functions

# ggplot2 theme for plotting
	theme_bw2 <- theme_set(theme_bw(base_size = 16))

# Function for calculating prediction intervals for plotting concentrations based on reactive input
# 90% prediction intervals
	CI80lo <- function(x) quantile(x, probs = 0.1)
	CI80hi <- function(x) quantile(x, probs = 0.9)
# 95% prediction intervals
	CI90lo <- function(x) quantile(x, probs = 0.05)
	CI90hi <- function(x) quantile(x, probs = 0.95)

# Time range - times where a concentration will be calculated
	set.time <- seq(from = 0, to = 120, by = 0.25)

# Define parameter values
# Thetas
	CLPOP <- 10  # Clearance, L/h
	V1POP <- 50  # Volume of central compartment, L
	QPOP <-  10  # Inter-compartmental clearance, L/h
	V2POP <- 100  # Volume of peripheral compartment, L
	KAPOP <- 0.5  # Absorption rate constant, h^-1

	COV1 <- 0.5  # Effect of smoking status
	COV2 <- 1.15  # Effect of creatinine clearance on clearance

# Omegas (as SD)
	ETA1SD <- 0.16  # PPV for clearance
	ETA2SD <- 0.16  # PPV for volume of central compartment
	ETA3SD <- 0.16  # PPV for absorption rate constant

	SDVAL <- c(ETA1SD, ETA2SD, ETA3SD)

# Specify a correlation matrix for ETA's
	R12 <- 0.5  # Correlation coefficient for CL-V1
	R13 <- 0.7  # Correlation coefficient for CL-KA
	R23 <- 0.5  # Correlation coefficient for V1-KA

# Calculate ETA values for each subject
	cor.vec <- c(
		1, R12, R13,
		R12, 1, R23,
		R13, R23, 1)
	CORR <- matrix(cor.vec, 3, 3)

# Use this function to turn CORR and SDVAL into a covariance matrix
	OMEGA <- cor2cov(CORR, SDVAL)

# Define continuous infusion
# This uses the approxfun() function
# Interpolates infusion rate for use in the differential equations
	inf.rate.fun <- function(time, rate) {
		approxfun(time, rate, method = "const")
	}

# Function containing differential equations for amount in each compartment
	DES <- function(T, A, THETA, inf.rate.fun) {
	# Determine infusion rate for value of T (time)
	  RateC <- inf.rate.fun(T)

	# Determine micro rate constants
	  K12 <- THETA[1]
	  K21 <- THETA[2]
	  K10 <- THETA[3]
		KA <- THETA[4]

	# Rate of change of A (amount) with respect to T (time)
	  dAdT <- vector(length = 3)
	  	dAdT[1] =       - KA*A[1]  # Depot - dose enters the system here
	  	dAdT[2] = RateC + KA*A[1] - K12*A[2] + K21*A[3] - K10*A[2]  # Central
	  	dAdT[3] =                   K12*A[2] - K21*A[3]  # Peripheral

	  list(dAdT)
	}

# Compile DES function
# It's called by lsoda for each individual in the dataset
  DES.cmpf <- cmpfun(DES)

# Function for simulating concentrations for each parameter set
  simulate.conc <- function(par.data, event.data, times, inf.rate.fun) {

  # List of parameters from input for the differential equation solver
    theta.list <- c(
			"K12" = par.data$K12,
			"K21" = par.data$K21,
			"K10" = par.data$K10,
			"KA" = par.data$KA)

  # Set initial compartment conditions
    A_0 <- c(A1 = 0, A2 = 0, A3 = 0)

  # Run differential equation solver for simulated variability data
    var.data <- lsoda(y = A_0, times = times, func = DES.cmpf, parms = theta.list,
			events = list(data = event.data), inf.rate.fun = inf.rate.fun)
    var.data <- as.data.frame(var.data)
  }

# Compile simulate.conc function
# It's called by ddply for each individual in the dataset
  simulate.conc.cmpf <- cmpfun(simulate.conc)
