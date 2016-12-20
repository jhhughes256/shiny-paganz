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
			h4("Patient Information"),
		# Slider input for age
			numericInput("age",
				"Age (years):",
				min = 0, max = 100, value = 32
			),  #numericInput
		# Slider input for weight
			numericInput("wt",
				"Total body weight (kg):",
				min = 0, max = 150, value = 70
			),  #numericInput
		# Select input for gender

		# Numeric input for serum creatinine

		# Display creatinine clearance from input

			br()

		# Radio buttons for smoking status

    # Heading for simulation options and horizontal line

		# Radio buttons for number of individuals

		# Checkbox for logscale

		# Select input for prediction intervals

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
				# IV Infusion Dosing
					checkboxInput("infcheck", "IV Infusion Dosing", FALSE),
					conditionalPanel(condition = "input.infcheck == true",
					# Slider input for IV infusion dose
						sliderInput("infdose",
							"Dose (mg):",
							min = 0, max = 1000, value = 500, step = 100
						),  #sliderInput
					# Slider input for IV infusion duration
						radioButtons("infdur",
							"Duration (hours):",
							choices = list(
								"2 hours" = 2,
								"12 hours" = 12),
							selected = 2,
							inline = TRUE
						),  #sliderInput
					# Numeric input for IV infusion starting time
						numericInput("inftimes",
							"Start Time (hours):",
							min = 0, max = 120, value = 72, step = 1
						),  #numericInput
					# Text output for IV infusion rate
						textOutput("infrate")
					)	#conditionalPanel
				)	#column

				# Oral dosing column

				# Checkbox for oral dosing

				# Conditional panel if checkbox is selected for oral dosing

				# Slider input for oral dose

				# Select input for oral dose frequency

				# Numeric input for start time of oral dosing

				### Oral end ###

				# IV dosing column

				# Checkbox for IV dosing

				# Conditional panel if checkbox is selected for IV dosing

				# Slider input for IV bolus dose

				# Slider input for IV bolus dose time

				### IV end ###

			)	#fixedRow
		)	 #mainPanel
	)	#sidebarLayout
)	#fixedPage
