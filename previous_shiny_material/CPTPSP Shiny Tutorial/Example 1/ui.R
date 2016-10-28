#Define UI for Flip-Flop PK Application (Example 1)
fluidPage(
	fluidRow(
	
		#Application Title
		h2("Flip-Flop Kinetics"),
		
		hr(),	#Add a break with a horizontal line
		
		#Plot output	
		plotOutput("plotCONC"),
    
		hr(),
    
		#Application Description
		p("For a one compartment first-order absorption model, there are two parameter sets that give identical results."),
		p("When elimination is rate limiting (KA = 2, KE = 0.2, V = 10) and when absorption is rate limiting (KA = 0.2, KE = 2, V = 1)."),
		h5("Show this by moving the sliders:"),
		
		br(),	#Add a blank break between widgets
		
		#Select values for rate constants using sliders
		sliderInput("KA",
					"Absorption Rate Constant (KA):",
					min = 0.2,
					max = 2,
					value = 2,
					step = 0.1),
					
		br(),
					
		sliderInput("KE",
					"Elimination Rate Constant (KE):",
					min = 0.2,
					max = 2,
					value = 0.2,
					step = 0.1),
					
		br(),			
					
		sliderInput("V",
					"Volume of Distribution (V):",
					min = 1,
					max = 10,
					value = 10,
					step = 0.05),			
					
	align = "center")	#Brackets closing "fluidRow"			
		
)	#Brackets closing "fluidPage"
