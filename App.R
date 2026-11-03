library(shiny)
library(shinyFiles)
library(htmltools)
library(bslib)
library(duckdb)
library(dplyr)
library(DBI)
library(DT)

# Source all R scripts in the R/ directory
invisible(lapply(list.files("src", pattern = "\\.R$", full.names = TRUE), source))

ui <- page_navbar(
  id = "main_nav",
  title = "MitoVarxploreR",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel("Setup / Load Data", value = "setup", setupUI("setup")),
  nav_panel("All Variants", value = "all_variants", allVariantsUI("all_variants")),
  nav_panel("Common Variants", value = "common_variants", commonVariantsUI("common_variants"))
)

server <- function(input, output, session) {
  # Initialize a single DuckDB connection for the entire session
  con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  # Ensure connection is closed when user leaves
  session$onSessionEnded(function() {
    dbDisconnect(con, shutdown = TRUE)
  })
  
  # 1. Call Setup Module. Returns a reactive containing the Parquet directory path.
  # We pass 'session' so the module can trigger a tab change in the main_nav.
  db_path <- setupServer("setup", parent_session = session)
  
  # 2. Call Data Modules, passing the reactive path and DB connection
  allVariantsServer("all_variants", db_path, con)
  commonVariantsServer("common_variants", db_path, con)
}

shinyApp(ui, server)