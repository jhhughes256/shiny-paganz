#Define UI for Preliminary Patient Education Tool Application (Example 2)
fixedPage(

	#Logo and Application Title
	fixedRow(
		column(10,
		h2("How many times a day and at what dose do you need to take this
		medication so that it is effective but not toxic?", align = "center"), offset = 1)	
		),	#Brackets closing "fixedRow"
				
	hr(),	#Add a break with a horizontal line
	
	#Sidebar Panel with Widgets
	sidebarLayout(
		sidebarPanel(
		
			#Heading
			h4("Patient Information"),
			
			#Slider input for patient weight
			sliderInput("WT",
						"Weight (kg):",
						min = 50,
						max = 120,
						value = 70,
						step = 5),
						
			#Slider input for patient age
			sliderInput("AGE",
						"Age (years):",
						min = 30,
						max = 90,
						value = 60,
						step = 30),

			br(), #Add a blank space between sections
			
			#Heading
			h4("Dosing Information"),
						
			#Slider input for prescribed dose
			sliderInput("DOSE",
						"Prescribed dose (mg):",
						min = 0,
						max = 1000,
						value = 500,
						step = 50),
      
      br(),
		
			#Selection box for dosing regimen
			selectInput("FREQ",
						"Dose Frequency:",
						choices = list("Once daily" = 1,
										"Twice daily" = 2,
										"Three times daily" = 3),
						selected = 1),
				
			#Question for missing or doubling doses
			h5("If you miss a dose, can you double the next one?", align = "left"),
												
				conditionalPanel(condition = "input.FREQ == 1",
				
				#Checkbox input for missing/doubling doses
				#Once daily dosing
				checkboxInput("DOSE13",
							"Dose Missed on Day 3:",
							value = FALSE),
							
				checkboxInput("DOSE14",
							"Dose Doubled on Day 4:",
							value = FALSE)
							
				),	#Brackets closing "conditionalPanel"
						
				conditionalPanel(condition = "input.FREQ == 2",
				
				#Checkbox input for missing/doubling doses
				#Twice daily dosing
				checkboxInput("DOSE26",
							"Second Dose Missed on Day 2:",
							value = FALSE),
							
				checkboxInput("DOSE27",
							"First Dose Doubled on Day 3:",
							value = FALSE)
				),	#Brackets closing "conditionalPanel"
				
				conditionalPanel(condition = "input.FREQ == 3",
				
				#Checkbox input for missing/doubling doses
				#Three times daily dosing
				checkboxInput("DOSE312",
							"Third Dose Missed on Day 3:",
							value = FALSE),
							
				checkboxInput("DOSE313",
							"First Dose Doubled on Day 4:",
							value = FALSE)					
							
				),	#Brackets closing "conditionalPanel"
			
			#Text
			h6("Helpful hint: refresh browser to reset values")
		
		),	#Brackets closing "siderbarPanel"
		
		#Plot output	
		mainPanel(
			plotOutput("plotCONC", height = 600, width = 800)
			
		)	#Brackets closing "mainPanel"
			
	)	#Brackets closing "sidebarLayout"
		
)	#Brackets closing "fixedPage"