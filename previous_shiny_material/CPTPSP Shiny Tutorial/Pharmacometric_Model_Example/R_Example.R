#Example R script for simulating a population and concentrations as described by:
#2-compartment model
#Dosing regimen: IV bolus (at 0.5 hours)
#				 IV infusion (at 0 hours, duration = 12 hours)
#				 IV bolus (at 15 hours)
#Population parameter variability on and CL, V1, V2 
#Variance-covariance matrix for  and CL, V1, V2
#Covariate effects: Weight (WT) on CL
#					Gender (SEX) on V1
#Proportional and additive residual error model


#---------------------------------------------------------------------------------------------------------------------------------------------
#Remove all current objects in the workspace
	rm(list=ls(all=TRUE))

#Load package libraries
	library(deSolve)	#Differential equation solver
	library(ggplot2)	#Plotting
	library(plyr)		#Split and rearrange data, ddply function
	library(grid)		#Plotting
	library(MASS)		#mvrnorm function
	library(MBESS)		#cor2cov function
	library(compiler)	#Compile repeatedly-called functions

	
#---------------------------------------------------------------------------------------------------------------------------------------------
#Step 1.  Create a parameter dataframe with ID and parameter values for each individual
	
#Define a population
	n <- 1000
	ID <- seq(from = 1, to = n, by = 1)
	WT <- rlnorm(n, meanlog = log(70), sdlog = 0.09)	#Log normal distribution for body weight (kg)
	SEX <- rbinom(n, size = 1, prob = 0.48)             #Binomial distribution for sex  0 is male, 1 is female  ; 48% females
	
	#Check demographics
	hist(WT)
	table(SEX)

	
#Define parameter values
    #Thetas	
   	CLPOP <- 10  	#Clearance L/hour
   	V1POP <- 50		#Central volume L
	QPOP <-  10  	#Inter-compartmental clearance L/h
	V2POP <- 100	#Peripheral volume L

	COV1 <- 0.75	#Covariate effect of body weight on clearance (power model)
	COV2 <- 0.5		#Covariate effect of female gender on V1

	#Omegas (as SD)
	ETA1SD <- 0.16
	ETA2SD <- 0.16
	ETA3SD <- 0.16
	
	#Specify a correlation matrix for ETA's
	R12 <- 0.5	#Correlation coefficient for CL-V1
	R13 <- 0.7	#Correlation coefficient for CL-V2
	R23 <- 0.5	#Correlation coefficient for V1-V2
	
	#Epsilons (as SD)
	EPS1SD <- 0.3  #Proportional residual error
    EPS2SD <- 0.1  #Additive residual error	

	
#Calculate ETA values for each subject
   	CORR <- matrix(c(1,R12,R13,R12,1,R23,R13,R23,1),3,3)
		
    #Specify the between subject variability for CL, V1, V2
	SDVAL <- c(ETA1SD,ETA2SD,ETA3SD)
		
    #Use this function to turn CORR and SDVAL into a covariance matrix
	OMEGA <- cor2cov(CORR,SDVAL)
		
    #Now use multivariate rnorm to turn the covariance matrix into ETA values
	ETAmat <- mvrnorm(n = n, mu = c(0,0,0), OMEGA)   
	ETA1 <- ETAmat[,1]
	ETA2 <- ETAmat[,2]
	ETA3 <- ETAmat[,3]


#Define individual parameter values
    CL <- CLPOP*exp(ETA1)*(WT/70)^COV1	
	V1 <- V1POP*exp(ETA2)*(1 + SEX*COV2)
	Q  <- QPOP
	V2 <- V2POP*exp(ETA3)
	
		
#Calculate rate-constants for differential equation solver
	K12 <- Q/V1
	K21 <- Q/V2
	K10 <- CL/V1

	
#Collect the individual parameter values in a parameter dataframe	
    par.data <- data.frame(ID,CL,V1,Q,V2,K12,K21,K10,WT,SEX)
	head(par.data)

	
#---------------------------------------------------------------------------------------------------------------------------------------------	
#Step 2.  Make a time sequence and specify the dose information for a system of differential equations
#There are a number of ways that doses can be coded for the deSolve package - see help for deSolve
	
#Creating continuous infusion (for 24 hours) - this use the approxfun function to make a "forcing function" for infusion rate in the differential equations
	CDOSE <- 650	#mg	
	CTinf <- 12 	#hours
	CRATE <- CDOSE/CTinf
			
	CTIMEinf <- c(0,CTinf,100) #100 - make sure the function works long after the 24 hours
	CRATEinf <- c(CRATE,0,0)
		
    #Define an interpolation function that returns rate when given time - "const"
    #i.e. at TIME = 0, RATE = CRATE. At TIME = 12, RATE = 0 (infusion stopped)
	Cstep.doseinf <- approxfun(CTIMEinf, CRATEinf, method = "const")
		#Testing
		Cstep.doseinf(0)   #Returns the infusion rate for time = 0 etc.
		Cstep.doseinf(11.9)
		Cstep.doseinf(12)
		Cstep.doseinf(12.1)
	
	
#Specify bolus doses - this uses the option to specify "events" in deSolve using a dataframe - see help for deSolve
	B1DOSE <- 250	#mg for first bolus dose
	B2DOSE <- 150	#mg for second bolus dose
	B1T <- 0.5        #time of first bolus (h)
	B2T <- 15       #time of second bolus (h)
	
    #Define bolus dose events							
	BOLUSDOSEdata <- data.frame(var = c(1,1),  #adding to "compartment" 1
				time = c(B1T,B2T),
				value = c(B1DOSE,B2DOSE),
				method = c("add","add"))
		
				
#Make a TIME sequence (hours)
	TIME <- seq(from = 0, to = 24, by = 0.25)
    TIME <- sort(unique(c(TIME,CTinf,B1T,B2T)))	 
	   #The time sequence must include all "event" times for deSolve so add them here and sort
	   #Do not repeat a time so use "unique" as well
			

#Function containing differential equations for amount in each compartment
	DES <- function(T, A, THETA) {
	
		RateC <- Cstep.doseinf(T)	#Infusion rate
		
		K12 <- THETA[1]                                   
		K21 <- THETA[2]
		K10 <- THETA[3]
			
		dA <- vector(length = 2)
        dA[1] =  RateC -K12*A[1] +K21*A[2] -K10*A[1]	#Central compartment 
        dA[2] =  K12*A[1] -K21*A[2]						#Peripheral compartment 
			
		list(dA)
	}

#Compile DES function - it's called by lsoda for each individual in the dataset	
	DES.cmpf <- cmpfun(DES)

			
#---------------------------------------------------------------------------------------------------------------------------------------------
#Step 3.  Run the differential equation solver for each patient in par.data
#Function for simulating concentrations for the ith patient
	simulate.conc <- function(par.data) {	
	
	#List of parameters from input for the differential equation solver			
		THETAlist <- c("K12" = par.data$K12,"K21 "= par.data$K21,"K10" = par.data$K10)
		
	#Set initial compartment conditions
		A_0 <- c(A1 = 0, A2 = 0)
		
	#Run differential equation solver for simulated variability data	
		var.data <- lsoda(A_0, TIME, DES.cmpf, THETAlist, events = list(data = BOLUSDOSEdata))
		var.data <- as.data.frame(var.data)	
	}

#Compile simulate.conc function	- it's called by ddply for each individual in the dataset	
	simulate.conc.cmpf <- cmpfun(simulate.conc)
	
#Apply simulate.conc.cmpf to each individual in par.data
#Whilst maintaining their individual values for V1, SEX and WT for later calculations
	sim.data <- ddply(par.data, .(ID,V1,SEX,WT), simulate.conc.cmpf)
	
#Calculate individual concentrations in the central compartment	
	sim.data$IPRE <- sim.data$A1/sim.data$V1

	
#---------------------------------------------------------------------------------------------------------------------------------------------
#Step 4.  Add residual error 
	
#Use random number generator to simulate residuals from a normal distribution
    nobs <- n*length(TIME)	 #number of observations = number of subjects * number of time points per subject
	EPS1 <- rnorm(nobs, mean = 0, sd = EPS1SD)	#Proportional residual error
	EPS2 <- rnorm(nobs, mean = 0, sd = EPS2SD)	#Additive residual error
	sim.data$DV <- sim.data$IPRE*(1 + EPS1) + EPS2
	
	
#---------------------------------------------------------------------------------------------------------------------------------------------
#Step 5.  Draw some plots of the simulated data

#Use custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 16))  
	theme_bw2 <- theme_update(plot.margin = unit(c(1.1,1.1,3,1.1), "lines"),
	axis.title.x=element_text(size = 16, vjust = 0),
	axis.title.y=element_text(size = 16, vjust = 0, angle = 90),
	strip.text.x=element_text(size = 14),
	strip.text.y=element_text(size = 14, angle = 90))

#Factor covariate values (for example, SEX)
	sim.data$SEXf <- as.factor(sim.data$SEX)
	levels(sim.data$SEXf) <- c("Male","Female")
	
#Function for calculating 5th and 95th percentiles for plotting concentrations
	CIlow <- function(x) quantile(x, probs = 0.05)
	CIhi <- function(x) quantile(x, probs = 0.95)	
	
#Generate a plot of the sim.data
	plotobj <- NULL
	plotobj <- ggplot(data = sim.data)
	plotobj <- plotobj + stat_summary(aes(x = time, y = DV), fun.ymin = CIlow, fun.ymax = CIhi, geom = "ribbon", fill = "blue", alpha = 0.2)
	plotobj <- plotobj + stat_summary(aes(x = time, y = IPRE), fun.ymin = CIlow, fun.ymax = CIhi, geom = "ribbon", fill = "red", alpha = 0.3)
	plotobj <- plotobj + stat_summary(aes(x = time, y = IPRE), fun.y = median, geom = "line", size = 1, colour = "red")
	plotobj <- plotobj + scale_y_continuous("Concentration (mg/L) \n")
	plotobj <- plotobj + scale_x_continuous("\nTime (hours)")
	print(plotobj)
	
#Facet for SEX	
	plotobj1 <- plotobj + facet_wrap(~SEXf)
	print(plotobj1)
	
