# Define UI
fixedPage(
# Application Title and Logo
	fixedRow(
		column(10,
			h2("PAGANZ - Shiny Example App"),
			h6(""),
    	offset = 1, align = "center"
		)  #column
	),	#fixedRow

# Add a break with a horizontal line
	hr(),

# Sidebar panel with widgets
	sidebarLayout(
		sidebarPanel(
		# Oral Dosing
			checkboxInput("pocheck", "Oral Dosing", FALSE),
			conditionalPanel("input.pocheck == true",
			# Slider input for oral dose to be given every 24 hours
				sliderInput("podose",
					"Oral Dose (mg):",
					min = 0, max = 1000, value = 500, step = 1
				),  #sliderInput
			# Slider input for oral dose frequency
				sliderInput("potimes",
					"Frequency of Dosing (times per day):",
					min = 1, max = 4, value = 1, step = 1
				)  #sliderInput
			),  #conditionalPanel

		# IV Bolus Dosing
			checkboxInput("ivcheck", "IV Bolus Dosing"),
			conditionalPanel("input.ivcheck == true",
			# Slider input for IV bolus dose
				sliderInput("ivdose",
					"IV Bolus Dose (mg):",
					min = 0, max = 1000, value = 500, step = 1
				),  #sliderInput
			# Slider input for IV bolus dose time
				numericInput("ivtimes",
					"IV Bolus Dose Time (hours):",
					min = 0, max = 120, value = 0
				)  #numericInput
			),  #conditionalPanel

		# Infusion Dosing
			checkboxInput("infcheck", "Infusion Dosing"),
			conditionalPanel("input.infcheck == true",
			# Slider input for IV infusion dose
				sliderInput("infdose",
					"Infusion Dose (mg):",
					min = 0, max = 4000, value = 2000, step = 1
				),  #sliderInput
			# Slider input for IV infusion duration
				sliderInput("infdur",
					"Infusion Duration (hours):",
					min = 0, max = 24, value = 10, step = 1
				),  #sliderInput
			# Slider input for IV infusion starting time
				sliderInput("inftimes",
					"Infusion Start Time (hours):",
					min = 0, max = 120, value = 72, step = 1
				)  #sliderInput
			),  #conditionalPanel

		# Covariates
			checkboxInput("covcheck", "Covariate Information"),
			conditionalPanel("input.covcheck == true",
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
					"Serum creatinine (umol/L):",
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
				)  #radioButtons
			),  #conditionalPanel

		# Slider input for number of individuals
			sliderInput("n",
				"Number of Individuals:",
				min = 10, max = 1000, value = 10, step = 10
			),  #sliderInput

			br(),

		# Button to initiate simulation
			submitButton("Simulate"),
			align = "left"	 #sidebarPanel options
		),  #sidebarPanel

	# Main panel to contain the concentration versus time plot
		mainPanel(
		# Plot output for concentration-time profile
			plotOutput("plotconc", height = 650, width = 750),
			align = "center"  #mainPanel options
		)	 #mainPanel
	)	#sidebarLayout
)	#fixedPage
