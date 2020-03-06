
# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Bight Data"),
  sidebarLayout(sidebarPanel(
    selectInput("reference_metal",
                label = "Select a Reference metal",
                choices = calibration.env$reference_metals,
                selected = calibration.env$reference_metals[2],
    ),
    selectInput("trace_metal",
                label = "Select a Trace Metal",
                choices = calibration.env$trace_metals
    ),
    checkboxInput("debug_mode","Calibration Values"),
    actionButton("rm_button", "Reference Metal Overview"),
    actionButton("tm_button", "Trace Metal Detail"),
  ),
  mainPanel( 
    #this will create a space for us to display our map
    
    tabsetPanel(id="tabs",type = "tabs",
                tabPanel("Map", 
                         fluidRow(
                           leafletOutput(outputId = "mymap"),
                         ),
                         fluidRow(verbatimTextOutput("Click_text")),
                         ), 
                tabPanel("Calibration Curve", 
                         h3("Calibration Plot"),
                         plotOutput("calibrationPlot"),
                         h3("Residual Analysis"),
                         h4("Residuals By Fit"),
                         plotOutput("residualsPlot"),
                         h4("Residuals By Latitude"),
                         plotOutput("residualsByLat"),
                         h4("Residuals By Longitude"),
                         plotOutput("residualsByLong"),
                         h4("Ratio Actual to Predicted By Latitude"),
                         plotOutput("APByLat"),
                         h4("Ratio Actual to Predicted By Longitude"),
                         plotOutput("APByLong")
                         ),
                tabPanel("Table", h3("Model Summary"),
                         fluidRow(verbatimTextOutput("TraceModelSummary")),
                         h3("Predicted Values"),
                         fluidRow(tableOutput("TraceMetalPredictions")),
                ),
                tabPanel("Overview Table",h3("Normalized Summary"),
                         fluidRow(tableOutput("OverviewTableNormal")),
                         h3("Raw Summary"),
                         fluidRow(tableOutput("OverviewTableRaw"))
                )
    ),

  )
  ))