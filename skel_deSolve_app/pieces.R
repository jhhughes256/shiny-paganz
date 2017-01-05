# Covariate Server Code

  ##########
  ##_CRCL_##
  ##########

  # Add crcl reactive text output
  crcl <- ((140 - input$age)*input$wt)/(input$secr*0.815)
  if (input$sex == 0) crcl <- ((140 - input$age)*input$wt)/(input$secr*0.815)*0.85
  paste0("Creatinine clearance = ", signif(crcl, digits = 3), " mL/min")

# ------------------------------------------------------------------------------
# Simulation Options Server Code

  ########
  ##_CI_##
  ########

  # Code for switching between levels of prediction intervals
  if (input$ci == 1) {
    plotobj <- plotobj + stat_summary(aes(x = time, y = IPRED),
      fun.ymin = CI80lo, fun.ymax = CI80hi,
      geom = "ribbon", fill = "red", alpha = 0.3)
  } else {
    plotobj <- plotobj + stat_summary(aes(x = time, y = IPRED),
      fun.ymin = CI90lo, fun.ymax = CI90hi,
      geom = "ribbon", fill = "red", alpha = 0.3)
  }

  #########
  ##_LOG_##
  #########

  # Code for adding log scale to y axis
  if (input$logscale == FALSE) {
    plotobj <- plotobj + scale_y_continuous("Concentration (mg/L) \n",
      breaks = seq(from = 0, to = max(Rsim.data()$IPRED), by = ceiling(max(Rsim.data()$IPRED)/10)),
      lim = c(0, max(Rsim.data()$IPRED)))
  } else {
    plotobj <- plotobj + scale_y_log10("Concentration (mg/L) \n",
      breaks = c(0.01,0.1,1,10), labels = c(0.01,0.1,1,10),
      lim = c(NA, max(Rsim.data()$IPRED)))
  }

# ------------------------------------------------------------------------------
# Dosing Server Code

  ##########
  ##_ORAL_##
  ##########

# Specify oral doses
# This uses the option to specify "events" in deSolve using a dataframe
	if (input$pocheck == FALSE) {
		oral.dose <- 0
		oral.dose.times <- 0
	} else {
		oral.dose <- input$podose
		if (input$potimes == 1) pofreq <- 24
		if (input$potimes == 2) pofreq <- 12
		if (input$potimes == 3) pofreq <- 8
		if (input$potimes == 4) pofreq <- 6
		oral.dose.times <- seq(from = input$postart, to = 120, by = pofreq)
	}
# Define bolus dose events
# Below works for constant dosing
	oral.dose.data <- data.frame(
		var = 1,  # Enters into depot compartment (A[1])
    time = oral.dose.times,
    value = oral.dose,
    method = "add"
	)

  ########
  ##_IV_##
  ########

# Specify bolus intravenous doses
# Specifies "events" as seen in oral dosing
	if (input$ivcheck == FALSE) {
		iv.dose <- 0
		iv.dose.times <- 0
	} else {
		iv.dose <- input$ivdose  # First bolus dose (mg)
	  iv.dose.times <- input$ivtimes  # Time of first bolus (h)
	}

# Define bolus dose events
	iv.dose.data <- data.frame(
		var = 2,  # Enters into central compartment (A[2])
    time = iv.dose.times,
    value = iv.dose,
    method = "add")

# EXTRA

# Combine dose data
  all.dose.data <- rbind(oral.dose.data, iv.dose.data)
