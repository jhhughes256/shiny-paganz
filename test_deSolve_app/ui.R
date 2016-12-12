# Define UI
fixedPage(
# Application Title and Logo
	fixedRow(
		h2("PAGANZ - Shiny Example App"),
		align = "center"
	),	#fixedRow

# Add a break with a horizontal line
	hr(),

# Sidebar panel with widgets
	sidebarLayout(
		sidebarPanel(
		# Covariates
		# Slider input for age
			numericInput("age",
				"Age (years):",
				min = 0, max = 100, value = 32
			),  #numericInput
		# Slider input for weight
			numericInput("wt",
				"Total body weight (kg):",
				min = 0, max = 72, value = 12
			),  #numericInput
		# Slider input for serum creatinine
			numericInput("secr",
				"Serum creatinine (µmol/L):",
				min = 0, max = 72, value = 12, step = 6
			),  #numericInput
		# Radio buttons for gender
			selectInput("sex",
				"Gender:",
				choices = list(
					"Male" = 0,
					"Female" = 1),
				selected = 1
			),  #radioButtons
		# Radio buttons for smoking status
			radioButtons("smok",
				"Smoking Status:",
				choices = list(
					"Not a Current Smoker" = 0,
					"Current Smoker" = 1),
				selected = 1
			),  #radioButtons
		# Slider input for number of individuals
			sliderInput("n",
				"Number of Individuals:",
				min = 10, max = 1000, value = 10, step = 10
			)  #sliderInput
		),  #sidebarPanel

	# Main panel to contain the concentration versus time plot
		mainPanel(
		# Plot output for concentration-time profile
			fixedRow(
				plotOutput("plotconc")
			),	#fixedRow
			hr(),
			fixedRow(
				column(4,
				# Oral Dosing
					checkboxInput("pocheck", "Oral Dosing", value = FALSE),
					conditionalPanel(condition = "input.pocheck",
					# Slider input for oral dose to be given every 24 hours
						sliderInput("podose",
							"Dose (mg):",
							min = 0, max = 1000, value = 500, step = 100
						),  #sliderInput
					# Slider input for oral dose frequency
						selectInput("potimes",
							"Dosing Frequency:",
							choices = list(
								"Once daily (every 24 hours)" = 1,
								"Twice daily (every 12 hours)" = 2,
								"Three times a day (every 8 hours)" = 3,
								"Four times a day (every 6 hours)" = 4
							),
							selected = 1
						)  #selectInput
					)	#conditionalPanel
				),	#column
				column(4,
				# IV Bolus Dosing
					checkboxInput("ivcheck", "IV Bolus Dosing", FALSE),
					conditionalPanel(condition = "input.ivcheck == true",
					# Slider input for IV bolus dose
						sliderInput("ivdose",
							"Dose (mg):",
							min = 0, max = 1000, value = 500, step = 100
						),  #sliderInput
					# Slider input for IV bolus dose time
						numericInput("ivtimes",
							"Dose Time (hours):",
							min = 0, max = 120, value = 0
						)  #numericInput
					)	#conditionalPanel
				),	#column
				column(4,
				# IV Infusion Dosing
					checkboxInput("infcheck", "IV Infusion Dosing", FALSE),
					conditionalPanel(condition = "input.infcheck == true",
					# Slider input for IV infusion dose
						sliderInput("infdose",
							"Dose (mg):",
							min = 0, max = 4000, value = 2000, step = 400
						),  #sliderInput
					# Slider input for IV infusion duration
						radioButtons("infdur",
							"Duration (hours):",
							choices = list(
								"30 minutes" = 1,
								"2 hours" = 2),
							selected = 1,
							inline = TRUE
						),  #sliderInput
					# Slider input for IV infusion starting time
						numericInput("inftimes",
							"Start Time (hours):",
							min = 0, max = 120, value = 72, step = 1
						)  #sliderInput
					)	#conditionalPanel
				)	#column
			)	#fixedRow
		)	 #mainPanel
	)	#sidebarLayout
)	#fixedPage
