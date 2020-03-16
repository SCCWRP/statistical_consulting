
server <- function(input, output, session) {
  
  # set global data used by the app
  calibration_sites=F
  clean_sites = calibration.env$clean_sites
  clean_sites_normalized = calibration.env$clean_sites_normalized
  dirty_sites = calibration.env$dirty_sites
  #create the map
  output$mymap <- renderLeaflet({
    leaflet(data) %>% 
      setView(lng = -118.16, lat = 33.75, zoom = 7) #setting the view over ~ center of bight
  })
  

    
  # Reference Metals Overview pageset
  observeEvent(input$rm_button,{
    hideTab(inputId="tabs",target="Table")
    hideTab(inputId="tabs",target="Calibration Curve")
    showTab(inputId="tabs",target="Overview Table")
    reference=input$reference_metal
    
    calibration_sites = input$debug_mode

    if(input$debug_mode){
      site_used = clean_sites
    } else {
      site_used = dirty_sites
    }
    pal <- colorNumeric(
      palette = c("blue"),
      na.color = "red",
      domain = clean_sites[["PPH"]])
    
    proxy <- leafletProxy("mymap")
    
    proxy %>% setView(lng = -118.16, lat = 33.75, zoom = 7)  %>% #setting the view over ~ center of bight
      clearMarkers() %>%
      addTiles() %>% 
      addCircleMarkers(data = site_used, lat = ~ lat, lng = ~ long, layerId = ~stationid, radius=5,  fillOpacity = 0.5,color=~pal(site_used[["PPH"]])) %>%
      addProviderTiles(providers$CartoDB.Positron)
   

    
    
    
    #generate overview_tables
    
    kabel_raw = calibration.env$rm_model_summary_kable(reference,clean_sites,"ReferenceMetal","PPH","TraceMetal","PPM")
    kabel_normalized = calibration.env$rm_model_summary_kable(reference,clean_sites_normalized,"ReferenceMetal","PPH","TraceMetal","PPM")
    
    output$OverviewTableNormal = function(){
      kabel_normalized
    }
    output$OverviewTableRaw = function(){
      kabel_raw
    }
     
  })

  
  # Trace Metal Overview Pageset
  observeEvent(input$tm_button,{
    showTab(inputId="tabs",target="Table")
    hideTab(inputId="tabs",target="Overview Table")
    showTab(inputId="tabs",target="Calibration Curve")
    calibration_sites = input$debug_mode
    
    tryCatch({calibration_curve = calibration.env$calibration_plot(input$reference_metal,input$trace_metal,clean_sites_normalized,dirty_sites,calibration_sites)},
             error = function(e){
               calibration_curve <<- calibration.env$calibration_plot(input$reference_metal,input$trace_metal,clean_sites,dirty_sites,calibration_sites)
               calibration_curve[["pointsPlot"]]<<- calibration_curve[["pointsPlot"]]+labs(caption="NOTE: Could not normalize the baseline relationship. Prediction intervals may be skewed")
               
             })

    output$calibrationPlot = renderPlot({
      
      plot(calibration_curve[["pointsPlot"]])})
    
    output$residualsPlot = renderPlot({
      
      plot(calibration_curve[["residualsPlot"]])})
    
    
    output$residualsByLat = renderPlot({
      
      plot(calibration_curve[["residualsByLat"]])})
    
    output$residualsByDepth = renderPlot({
      
      plot(calibration_curve[["residualsByDepth"]])})
    
    
    output$residualsByLong = renderPlot({
      
      plot(calibration_curve[["residualsByLong"]])})
    
    output$APByLat = renderPlot({plot(calibration_curve[["APByLat"]])})
    output$APByLong = renderPlot({plot(calibration_curve[["APByLong"]])})
    
    
    reference=input$reference_metal
    
    if(input$debug_mode){
      site_used = clean_sites
    } else {
      site_used = dirty_sites
    }
    
    # print out the summary tables
    
    output$TraceModelSummary = renderPrint({calibration_curve[["model"]]})
    
    output$TraceMetalPredictions = function(){
      predict_vals = predict(calibration_curve[["model"]],newdata = calibration_curve[["dirty_sites"]])
      df = as.data.frame(calibration_curve[["dirty_sites"]]$PPM)
      colnames(df) = c(input$trace_metal)
      df[["Predicted Amount"]] = predict_vals
      df[["Station"]] = calibration_curve[["dirty_sites"]]$stationid
      
      df[["Stratum"]] = calibration_curve[["dirty_sites"]]$Stratum
      df = df[,c("Station","Stratum",input$trace_metal,"Predicted Amount")]
      kabel = kable(df, format ="html",booktabs = T,escape=F,align="c") %>%
             kable_styling(position = "center")
      return(kabel)
    }
    
    
    
    
    pal <- colorNumeric(
      palette = c("blue"),
      na.color = "red",
      domain = calibration_curve$predicted_PPM[["fit"]])
    
    proxy <- leafletProxy("mymap")

    proxy %>% setView(lng = -118.16, lat = 33.75, zoom = 7)  %>% #setting the view over ~ center of bight
      clearMarkers() %>%
      addTiles() %>% 
      addCircleMarkers(data = calibration_curve[["predicted_PPM"]], lat = ~ lat, lng = ~ long, radius=5, layerId = ~stationid,   fillOpacity = 0.5,color=~pal(site_used[["Actual"]])) %>%
      addProviderTiles(providers$CartoDB.Positron)
    
  })
  
  
  observe({
    proxy <- leafletProxy("mymap")
    click<-input$mymap_marker_click
    if(is.null(click))
      return()
    text<-paste(click$id)

    text2 <- "tableInfoHere"
    proxy%>%
      clearPopups()%>%
    addPopups(dat=click,lat= click$lat, lng= click$lng, text)
    
    output$Click_text<-renderText({
      text2
    })
  })
  
  
  
    
  
  
    
  # select trace metal and contaminant, spit out table with:
  # mse, residual, model summary, map of anthropogenic effect
  # table from the paper (slope, intercept), table 3
  
  
  
  
}
