# R script for simulating a population from a generic 2-compartment model
# ------------------------------------------------------------------------------
# Load package libraries
	library(ggplot2)	# Plotting
	library(grid)	# Plotting
	library(dplyr)	# New plyr - required for mrgsolve
	library(mrgsolve)	# Metrum differential equation solver for pharmacometrics
# Define a custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 16))

# ------------------------------------------------------------------------------
# Set number of individuals that make up the 95% prediction intervals
	n <- 1000
# 95% prediction interval functions - calculate the 2.5th and 97.5th percentiles
	CI95lo <- function(x) quantile(x,probs = 0.025)
	CI95hi <- function(x) quantile(x,probs = 0.975)
# 90% prediction interval functions - calculate the 5th and 95th percentiles
	CI90lo <- function(x) quantile(x,probs = 0.05)
	CI90hi <- function(x) quantile(x,probs = 0.95)
# Set seed for reproducible numbers
	set.seed(123456)

# Time
	time.function <- function(i) {
		TIME <- c(seq(from = 0,to = 3,by = 0.25),seq(from = 4,to = 24,by = 1))+i*24
	}
	TIME <- lapply(0:4, function(i) time.function(i))
	TIME <- unique(unlist(TIME))

	TIME <- seq(from = 0,to = 144,by = 0.25)

# ------------------------------------------------------------------------------
# Define the model parameters and equations
	# Using mrgsolve - analytical solutions
	# This compiled model is used for simulating n individuals and their concentration-time profiles
		code <- '
		$INIT			// Initial conditions for compartments
							DEPOT = 0,	// Depot - dose enters the system here
							CENT = 0,	// Central
							PERI = 0,	// Peripheral
							AUC = 0	// Area under the curve compartment

		$PARAM		// Population parameters
							POPCL = 10,	// Clearance, L/h
							POPV1 = 50,	// Volume of central compartment, L
							POPQ = 10,	// Inter-compartmental clearance, L/h
							POPV2 = 100,	// Volume of peripheral compartment, L
							POPKA = 0.5,	// Absorption rate constant, h^-1

							// Covariate effects
							COV1 = 0.5,	// Effect of current smoking on clearance
							COV2 = 1.15,	// Effect of creatinine clearance on clearance

							// Default covariate values for simulation
							WT = 70,	// Total body weight (kg)
							SEX = 0,	// Gender - Male (0), Female (1)
							AGE = 60,	// Age (years)
							SECR = 100, // Serum creatinine, micromol/L
							SMOK = 0	// Smoking status - Not current (0), Current (1)

		$OMEGA		block = TRUE
							labels = s(BSV_CL,BSV_V1,BSV_KA)
							0.0256	// BSV for CL
							0.0128 0.0256	// BSV for V1
							0.0179 0.0128 0.0256	// BSV for KA

		$SIGMA		block = FALSE
							labels = s(ERR_PRO)
							0.09	// Proportional error

		$MAIN			// Covariate effects
							double SMOKCOV = 1;	// Not current smoker
							if (SMOK == 1) SMOKCOV = 1 + COV1;	// Current smoker

							double CRCL = ((140-AGE)*WT)/(SECR*0.815);	// Male creatinine clearance
							if (SEX == 1) CRCL = ((140-AGE)*WT)/(SECR*0.815)*0.85;	// Female creatinine clearance

							// Individual parameter values
							double CL = POPCL*pow(WT/70,0.75)*pow(CRCL/90,COV2)*SMOKCOV*exp(BSV_CL);
							double V1 = POPV1*(WT/70)*exp(BSV_V1);
							double Q = POPQ*pow(WT/70,0.75);
							double V2 = POPV2*(WT/70);
							double KA = POPKA*exp(BSV_KA);

		$ODE			// Differential equations
							dxdt_DEPOT = -KA*DEPOT;
							dxdt_CENT = KA*DEPOT -Q/V1*CENT +Q/V2*PERI -CL/V1*CENT;
							dxdt_PERI = Q/V1*CENT -Q/V2*PERI;
							dxdt_AUC = CENT/V1;

		$TABLE		table(IPRE) = CENT/V1;
							table(DV) = table(IPRE)*(1+ERR_PRO);

		$CAPTURE	CL V1 Q V2 KA CRCL
		'
	# Compile the model code
		mod <- mcode("popDEMO",code)
		# There is opportunity to simply update model parameters after the model code has been compiled

# ------------------------------------------------------------------------------
# Simulate concentration-time profiles for the population
	ID <- 1:n
	ID2 <- sort(c(rep(ID,times = length(TIME))))
	time <- rep(TIME,times = length(ID))
	input.conc.data <- data.frame(
		ID = ID2,
		time,
		amt = 0,
		evid = 0,
		rate = 0,
		cmt = 1
	)

	oral.dose.times <- c(1,seq(from = 12,to = 120,by = 12))
	oral.dose.data <- input.conc.data[input.conc.data$time %in% oral.dose.times,]
	oral.dose.data$amt <- 500
	oral.dose.data$evid <- 1
	oral.dose.data$rate <- 0
	oral.dose.data$cmt <- 1

	iv.dose.times <- 0
	iv.dose.data <- input.conc.data[input.conc.data$time %in% iv.dose.times,]
	iv.dose.data$amt <- 500
	iv.dose.data$evid <- 1
	iv.dose.data$rate <- 0
	iv.dose.data$cmt <- 2

	inf.dose.times <- 72
	inf.dose.data <- input.conc.data[input.conc.data$time %in% inf.dose.times,]
	inf.dose.data$amt <- 2000
	inf.dose.data$evid <- 1
	inf.dose.data$rate <- 200
	inf.dose.data$cmt <- 2

	input.conc.data <- rbind(input.conc.data,oral.dose.data,iv.dose.data,inf.dose.data)
	input.conc.data <- input.conc.data[with(input.conc.data, order(input.conc.data$ID,input.conc.data$time)),]

	conc.data <- mod %>% data_set(input.conc.data) %>% mrgsim()
	conc.data <- as.data.frame(conc.data)	# Convert to a data frame so that it is more useful for me!

# ------------------------------------------------------------------------------
# Plot results
	plotobj1 <- NULL
	plotobj1 <- ggplot(conc.data)
	plotobj1 <- plotobj1 + stat_summary(aes(x = time,y = IPRE),geom = "line",fun.y = median,colour = "red",size = 1)
	plotobj1 <- plotobj1 + stat_summary(aes(x = time,y = IPRE),geom = "ribbon",fun.ymin = "CI90lo",fun.ymax = "CI90hi",fill = "red",alpha = 0.3)
	plotobj1 <- plotobj1 + scale_x_continuous("\nTime (hours)",lim = c(0,120))
	# plotobj1 <- plotobj1 + scale_y_log10("Concentration (mg/L)\n")
	plotobj1 <- plotobj1 + scale_y_continuous("Concentration (mg/L)\n",breaks = seq(from = 0,to = 30,by = 5),lim = c(0,25))
	print(plotobj1)
