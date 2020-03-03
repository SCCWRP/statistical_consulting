# anything that is required to be in the global namespace
source("global.R")

#user interface
source("ui.R")

# server/backend logic
source("server.R")

# Run the application 
shinyApp(ui = ui, server = server)
