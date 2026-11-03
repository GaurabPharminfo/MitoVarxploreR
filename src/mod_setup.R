library(shiny)
library(bslib)
library(shinyFiles)
library(fs) # Used for defining the home directory

setupUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$head(
      tags$style(HTML("
        /* Custom styling for the text outputs to make them look clean and code-like */
        .setup-output { font-family: monospace; color: #6c757d; padding-left: 15px; }
      "))
    ),
    
    h3("Data Setup", class = "mb-4"),
    
    card(
      card_header("Option A: Build New Database", class = "fw-bold"),
      card_body(
        # Directory containing sample folders
        div(class = "d-flex align-items-center mb-3",
            div(style = "width: 280px;", strong("Directory containing sample folders:")),
            shinyDirButton(ns("samples_dir_btn"), "Browse...", "Select Samples Directory", class = "btn btn-outline-primary btn-sm"),
            # Fixed textOutput by wrapping it in a span with the class
            span(class = "setup-output", textOutput(ns("samples_dir_out"), inline = TRUE))
        ),
        
        # Sample ID Mapping
        div(class = "d-flex align-items-center mb-3",
            div(style = "width: 280px;", strong("Sample ID Mapping:")),
            shinyFilesButton(ns("mapping_csv_btn"), "Browse...", "Select Sample Mapping", multiple = FALSE, accept = c(".csv"), class = "btn btn-outline-primary btn-sm"),
            span(class = "setup-output", textOutput(ns("mapping_csv_out"), inline = TRUE))
        ),
        
        # Output Directory for Parquet
        div(class = "d-flex align-items-center mb-4",
            div(style = "width: 280px;", strong("Output Directory for Parquet:")),
            shinyDirButton(ns("output_dir_btn"), "Browse...", "Select Output Directory", class = "btn btn-outline-primary btn-sm"),
            span(class = "setup-output", textOutput(ns("output_dir_out"), inline = TRUE))
        ),
        
        actionButton(ns("run_etl"), "Run ETL & Load Data", class = "btn-primary w-25")
      )
    ),
    
    card(
      card_header("Option B: Load Existing Database", class = "fw-bold"),
      card_body(
        # Existing Database Directory Chooser
        div(class = "d-flex align-items-center mb-4",
            div(style = "width: 280px;", strong("Path to existing Parquet Database:")),
            shinyDirButton(ns("existing_dir_btn"), "Browse...", "Select Existing Parquet Directory", class = "btn btn-outline-success btn-sm"),
            span(class = "setup-output", textOutput(ns("existing_dir_out"), inline = TRUE))
        ),
        
        actionButton(ns("load_existing"), "Load Existing Data", class = "btn-success w-25")
      )
    )
  )
}

setupServer <- function(id, parent_session) {
  moduleServer(id, function(input, output, session) {
    # Reactive value to store the finalized path to pass to other modules
    db_path <- reactiveVal(NULL)
    
    # Define root volumes for the file/folder choosers (Home directory + OS Drives)
    roots <- c(Home = fs::path_home(), shinyFiles::getVolumes()())
    
    # 1. Initialize the shinyFiles observers
    shinyDirChoose(input, "samples_dir_btn", roots = roots, session = session)
    shinyFileChoose(input, "mapping_csv_btn", roots = roots, session = session)
    shinyDirChoose(input, "output_dir_btn", roots = roots, session = session)
    shinyDirChoose(input, "existing_dir_btn", roots = roots, session = session)
    
    # 2. Parse paths safely (returns path string if selected, else NULL)
    samples_dir_path <- reactive({
      if (is.integer(input$samples_dir_btn)) return(NULL)
      parseDirPath(roots, input$samples_dir_btn)
    })
    
    mapping_csv_path <- reactive({
      if (is.integer(input$mapping_csv_btn)) return(NULL)
      parseFilePaths(roots, input$mapping_csv_btn)$datapath
    })
    
    output_dir_path <- reactive({
      if (is.integer(input$output_dir_btn)) return(NULL)
      parseDirPath(roots, input$output_dir_btn)
    })
    
    existing_dir_path <- reactive({
      if (is.integer(input$existing_dir_btn)) return(NULL)
      parseDirPath(roots, input$existing_dir_btn)
    })
    
    # 3. Render the selected paths in the UI
    output$samples_dir_out <- renderText({
      path <- samples_dir_path()
      if (is.null(path)) "No directory selected" else path
    })
    
    output$mapping_csv_out <- renderText({
      path <- mapping_csv_path()
      if (is.null(path) || length(path) == 0) "No file selected" else path
    })
    
    output$output_dir_out <- renderText({
      path <- output_dir_path()
      if (is.null(path)) "No directory selected" else path
    })
    
    output$existing_dir_out <- renderText({
      path <- existing_dir_path()
      if (is.null(path)) "No directory selected" else path
    })
    
    # 4. Action: Run ETL
    observeEvent(input$run_etl, {
      s_dir <- samples_dir_path()
      o_dir <- output_dir_path()
      map_path <- mapping_csv_path()
      
      # Ensure mandatory directories are selected
      if (is.null(s_dir) || is.null(o_dir)) {
        showNotification("Please select both a Samples Directory and an Output Directory.", type = "warning")
        return()
      }
      
      # Handle optional mapping CSV safely
      map_path_safe <- if (!is.null(map_path) && length(map_path) > 0) map_path else ""
      
      showNotification("Running ETL... This may take a while.", type = "message", duration = NULL, id = "etl_notif")
      
      success <- tryCatch({
        run_etl_pipeline(s_dir, map_path_safe, o_dir)
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        FALSE
      })
      
      removeNotification("etl_notif")
      if (success) {
        showNotification("ETL Complete! Loading data...", type = "message")
        db_path(o_dir)
        nav_select("main_nav", "all_variants", session = parent_session) # Route to next page
      }
    })
    
    # 5. Action: Load Existing Database
    observeEvent(input$load_existing, {
      e_dir <- existing_dir_path()
      
      if (is.null(e_dir)) {
        showNotification("Please select an existing database directory.", type = "warning")
        return()
      }
      
      if (dir.exists(e_dir)) {
        db_path(e_dir)
        nav_select("main_nav", "all_variants", session = parent_session) # Route to next page
      } else {
        showNotification("Directory not found on disk.", type = "error")
      }
    })
    
    return(db_path)
  })
}